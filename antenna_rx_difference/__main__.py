import itertools
import datetime
import struct
import pickle
import sys
import re
import enum
import tempfile
import time
import json
import os
import shelve
import glob
import math

from astropy.coordinates import EarthLocation
import astropy.units as u
from astroplan import Observer

import pyproj
import numpy as np
import matplotlib.pyplot as plt

import geopy
import geopy.distance

import maidenhead

GNUPLOT_COLORS = ['#00e673', '#e62e00', 'purple', '#4d4dff', '#2eb8b8']

MAIDENHEAD_REGEX = re.compile('^[A-Z]{2}[0-9]{2}$')

CALLSIGN_REGEX = re.compile('^[A-Z,0-9]{3,}')

class SignalObservation(object):
  __slots__ = ('ts', 'freq', 'direction', 'mode', 'snr', 'freq_offset', 'msg',)
  def __init__(self, ts, freq, direction, mode, snr, _, freq_offset, msg):
    self.ts = ts
    self.freq = freq
    self.direction = direction
    self.mode = mode
    self.snr = snr
    self.freq_offset = freq_offset
    self.msg = msg

  def __hash__(self):
    return hash((self.ts, self.freq,  self.msg,))   

  def get_tx_station(self):
    data = self.msg.split(' ')
    if len(data) < 3:
      return
    callsign = data[1]
    if CALLSIGN_REGEX.match(callsign) is None:
      return
    return callsign

  def get_maidenhead(self):
    data = self.msg.split(' ')
    if len(data) < 3: # a message should have at least pieces of information, if it includes maidenhead
      return
     
    last_data = data[-1]
    if last_data == 'RR73': # commonly send to indicate the communication is over
      return

    if MAIDENHEAD_REGEX.match(last_data):
      return last_data

class DataMode(enum.IntEnum):
  FT8 = 1

class SignalDirection(enum.IntEnum):
  RX = 1
  TX = 2

WSJTX_LINE_REGEX = '^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(.*)$'
def parse_wsjtx_output(fin):
  pattern = re.compile(WSJTX_LINE_REGEX)
  for line in fin:
    m = pattern.match(line)
    if m is None:
      continue

    ts, freq, direction, mode, snr, _, freq_offset, msg = m.groups()

    dot = freq.find('.')
    freq_scale = 1
    if dot != -1:
      dot_dist = len(freq) - dot
      freq_scale = 10**(dot_dist-1)
     

    freq = int(freq.replace('.', '')) * freq_scale
    freq_offset = int(freq_offset)
    snr = int(snr)
    
    mode = mode.lower()
    if mode == 'rx':
      mode = SignalDirection.RX
    elif mode == 'tx':
      mode = SignalDirection.TX
    msg = msg.strip()

    year = int(ts[0:2]) + 2000
    month = int(ts[2:4])
    day = int(ts[4:6])

    hour = int(ts[7:9])
    minute = int(ts[9:11])
    seconds = int(ts[11:])

    ts = datetime.datetime(year, month, day, hour, minute, seconds)
    yield SignalObservation(ts, freq, direction, mode, snr, None, freq_offset, msg)


tmpA = tempfile.TemporaryFile(mode='w+b')
tmpB = tempfile.TemporaryFile(mode='w+b')

def write_object_sequence(fin, objs):
  for x in objs:
    pickled = pickle.dumps(x)
    header = struct.pack('>H', len(pickled))
    fin.write(header)
    fin.write(pickled)

def read_object_sequence(fin):
  fin.seek(0)
  while True:
    header = fin.read(2)
    if len(header) == 0:
      break
    data_len, = struct.unpack('>H', header)
    yield pickle.loads(fin.read(data_len))

nameA = os.path.split(sys.argv[1])[-1].split('.')[0].replace('_', ' ')
nameB = os.path.split(sys.argv[2])[-1].split('.')[0].replace('_', ' ')

with open(sys.argv[1]) as fin:
  write_object_sequence(tmpA, parse_wsjtx_output(fin))

with open(sys.argv[2]) as fin:
  write_object_sequence(tmpB, parse_wsjtx_output(fin))
  
del fin
observer = Observer(location = EarthLocation.from_geodetic(0.0 * u.deg, 0.0 * u.deg, 0.0 * u.m))
#sunrise = observer.sun_rise_time()
#sunset = observer.sun_set_time()

def truncate_ts_to_hour_unix(ts):
  truncated_ts = datetime.datetime(ts.year, ts.month, ts.day, ts.hour)
  unix_ts = int(time.mktime(truncated_ts.timetuple()))
  return unix_ts

def find_time_range(inputA, inputB):
  earliest = None
  latest = None

  for obj in itertools.chain(*(read_object_sequence(x) for x in (inputA, inputB,))):
    if earliest is None or earliest < obj.ts:
      earliest = obj.ts

    if latest is None or latest > obj.ts:
      latest = obj.ts

  return (earliest, latest,)

earliest_time, latest_time = find_time_range(tmpA, tmpB)

def count_by_freq_and_hour(fin):
  data = {}
  for obj in read_object_sequence(fin):
    unix_ts = truncate_ts_to_hour_unix(obj.ts)
    k = (unix_ts, obj.freq,)
    v = data.get(k)
    if v is None:
      v = 0
    data[k] = v + 1
  return data

def find_all_hours_and_freqs(dataA, dataB):
  all_freqs = set()
  [all_freqs.add(x[1]) for x in dataA.keys()]
  [all_freqs.add(x[1]) for x in dataB.keys()]

  all_hours = set()
  [all_hours.add(x[0]) for x in dataA.keys()]
  [all_hours.add(x[0]) for x in dataB.keys()]
  all_hours = sorted(list(all_hours))

  freq_ordering = sorted(list(all_freqs))
  return (all_hours, freq_ordering,)

GNUPLOT_TIMEFMT_HRS = "set timefmt \"%Y/%m/%d %H\"\n"
GNUPLOT_X_FORMAT_HRS = "set format x \"%m-%d\\n%H\"\n"
DATETIME_FORMAT_FOR_GNUPLOT = '%Y/%m/%d %H'

def write_files_for_counts_by_hour(nameA, nameB, tmpA, tmpB):
  dataA = count_by_freq_and_hour(tmpA)
  dataB = count_by_freq_and_hour(tmpB)
  all_hours, all_freqs = find_all_hours_and_freqs(dataA, dataB)

  csv_filename = 'counts_by_hour.csv'
  with open(csv_filename, 'w') as fout:
    fout.write("# input file A: %s\n" % (nameA,))
    fout.write("# input file B: %s\n" % (nameB,))
    fout.write("# hours covered: %d\n" % (len(all_hours),) )
    fout.write("# frequencies covered: %d\n" % (len(all_freqs),) )
    fout.write("# ")
    cols = ['hour (unix ts)']
    
    for freq in all_freqs:
      cols.append('%d Hz count (%s)' % (freq, nameA,))
      cols.append('%d Hz count (%s)' % (freq, nameB,))
   
    fout.write(','.join(cols))

    fout.write("\n")

    for hour in all_hours:
      row = [datetime.datetime.utcfromtimestamp(hour).strftime(DATETIME_FORMAT_FOR_GNUPLOT)]
      for freq in all_freqs:
        k = (hour, freq,)
        v_a = dataA.get(k, 0)
        v_b = dataB.get(k, 0)
        row.append('%.1f' % (v_a,))
        row.append('%.1f' % (v_b,))

      fout.write(','.join(str(x) for x in row))
      fout.write("\n")

  with open('counts_by_hour.gnuplot', 'w') as fout:
    fout.write("set datafile separator comma\n")
    fout.write("set datafile missing \"NaN\"\n")

    for i, freq in enumerate(all_freqs):
      filename = "counts_by_hour_%dHz.png" % (freq,)
      freq_human = "%.3f MHz" % (freq/(10**6.0),)

      fout.write("set terminal pngcairo font \"Arial, 16\" size 1920, 1080\n")
      fout.write("set output \"%s\"\n" % (filename,))

      fout.write("set title \"Signals received per hour (%s)\"\n" % (freq_human,))
      fout.write("set xlabel \"Hour\"\n")
      fout.write("set ylabel \"Count\"\n")

      a_column_idx = 1 + (i*2)
      b_column_idx = 1 + (i*2 + 1)

      fout.write("set xdata time\n")
      fout.write(GNUPLOT_TIMEFMT_HRS)
      fout.write(GNUPLOT_X_FORMAT_HRS)
      fout.write("set grid\n")

      # TODO colors that aren't boring
      fout.write("plot \"%s\" using 1:%d with linespoints lw 3 pt 7 ps 2 lc rgb \"black\" title \"%s\", \\\n" % (csv_filename, a_column_idx + 1, nameA,))
      fout.write("     \"%s\" using 1:%d with linespoints lw 3 pt 7 ps 2 lc rgb \"red\" title \"%s\"\n" % (csv_filename, b_column_idx + 1, nameB,))
      fout.write("\n")


def maximum_distance_by_hour(rx_maidenhead, station_locator, fin):
  rx_coords = maidenhead.to_location(rx_maidenhead)
  data = {}
  for obj in read_object_sequence(fin):
    loc = obj.get_maidenhead()
    if loc is None:   
      callsign = obj.get_tx_station()
      if callsign is None:
        continue # no callsign, skip it
      loc = station_locator.lookup(callsign)
      if loc is None: # still no location, skip it
        continue 
    if len(loc) == 4:
      loc = loc + 'mm' # pick the middle

    tx_coords = maidenhead.to_location(loc)
    # TODO switch to pyproj
    dist_km = geopy.distance.great_circle(rx_coords, tx_coords).km
    unix_ts = truncate_ts_to_hour_unix(obj.ts)
    k = (unix_ts, obj.freq,)

    v = data.get(k)
    if v is None or v < dist_km:      
      data[k] = dist_km
  return data
    
def write_maximum_distance_by_hour(rx_maidenhead, station_locator, nameA, nameB, tmpA, tmpB):
  dataA = maximum_distance_by_hour(rx_maidenhead, station_locator, tmpA)
  dataB = maximum_distance_by_hour(rx_maidenhead, station_locator, tmpB)
  all_hours, all_freqs = find_all_hours_and_freqs(dataA, dataB)
   
  csv_filename = 'max_distance_by_hour.csv'
  with open(csv_filename, 'w') as fout:
    fout.write("# input file A: %s\n" % (nameA,))
    fout.write("# input file B: %s\n" % (nameB,))
    fout.write("# hours covered: %d\n" % (len(all_hours),) )
    fout.write("# frequencies covered: %d\n" % (len(all_freqs),) )
    fout.write("# ")
    cols = ['hour (unix ts)']
    
    for freq in all_freqs:
      cols.append('%d Hz max. km. (%s)' % (freq, nameA,))
      cols.append('%d Hz max. km. (%s)' % (freq, nameB,))
   
    fout.write(','.join(cols))

    fout.write("\n")

    for hour in all_hours:
      row = [datetime.datetime.utcfromtimestamp(hour).strftime(DATETIME_FORMAT_FOR_GNUPLOT)]
      for freq in all_freqs:
        k = (hour, freq,)
        v_a = dataA.get(k)
        v_b = dataB.get(k)
        if v_a is None:
          row.append('NaN')
        else:
          row.append('%.1f' % (v_a,))

        if v_b is None:
          row.append('NaN')
        else:
          row.append('%.1f' % (v_b,))

      fout.write(','.join(str(x) for x in row))
      fout.write("\n")

  with open('max_distance_by_hour.gnuplot', 'w') as fout:
    fout.write("set datafile separator comma\n")
    fout.write("set datafile missing \"NaN\"\n")
    for i, freq in enumerate(all_freqs):
      filename = "max_dist_by_hour_%dHz.png" % (freq,)
      freq_human = "%.3f MHz" % (freq/(10**6.0),)

      fout.write("set terminal pngcairo font \"Arial, 16\" size 1920, 1080\n")
      fout.write("set output \"%s\"\n" % (filename,))

      fout.write("set title \"Maximum reception range by hour (%s)\"\n" % (freq_human,))
      fout.write("set xlabel \"Hour\"\n")
      fout.write("set ylabel \"Distance (km.)\"\n")

      a_column_idx = 1 + (i*2)
      b_column_idx = 1 + (i*2 + 1)

      fout.write("set xdata time\n")
      fout.write(GNUPLOT_TIMEFMT_HRS)
      fout.write(GNUPLOT_X_FORMAT_HRS)
      fout.write("set grid\n")

      # TODO colors that aren't boring
      fout.write("plot \"%s\" using 1:%d with linespoints lw 3 pt 7 ps 2 lc rgb \"black\" title \"%s\", \\\n" % (csv_filename, a_column_idx + 1, nameA,))
      fout.write("     \"%s\" using 1:%d with linespoints lw 3 pt 7 ps 2 lc rgb \"red\" title \"%s\"\n" % (csv_filename, b_column_idx + 1, nameB,))
      fout.write("\n")

class StationLocator(object):
  def __init__(self): 
    tmpname = 'tmp.stationlocation'
    if len(tmpname) > 3:
       for fname in glob.glob(tmpname + '*'):
         print('removing ' + tmpname)
         os.unlink(fname)
    self.db = shelve.open(tmpname)

  def insert(self, callsign, gridsquare):
    self.db[callsign] = gridsquare

  def lookup(self, callsign):
    return self.db.get(callsign)

class SignalMatcher(object):
  def __init__(self):
     tmpname = 'tmp.signalmatcher'
     if len(tmpname) > 3:
       for fname in glob.glob(tmpname + '*'):
         print('removing ' + tmpname)
         os.unlink(fname)
     self.db = shelve.open(tmpname)

  def insert(self, cat, obj):
    k = '%x' % hash(obj)
    v = self.db.get(k)
    if v is None:
      v = {}

    sub_v = v.get(cat)
    if sub_v is None:
      sub_v = []
      v[cat] = sub_v
    sub_v.append(obj)
    self.db[k] = v

  def generate_matches(self):
    freq_tolerance = 100
    for k in self.db.keys():
      matched = set()

      data_for_k = self.db[k]
      
      categories = list(data_for_k.keys())
      for cat in categories:        
        entries = data_for_k[cat]
        
        for i, entry in enumerate(entries):
          # Search every other entry for matches
          for other_cat in categories:
            if other_cat == cat:
              continue
            for j, other_entry in enumerate(data_for_k[other_cat]):
              if entry.mode == other_entry.mode and other_entry.freq == other_entry.freq and other_entry.ts == other_entry.ts and entry.msg == other_entry.msg:
                freq_offset_delta = abs(entry.freq_offset - other_entry.freq_offset)
                if freq_offset_delta < freq_tolerance: # two receivers may not indicate the same freq offset exactly
                  match_k = tuple(sorted( [ (cat, i,), (other_cat, j,) ]))
                  if match_k not in matched:
                    matched.add(match_k)
                    match = { cat: entry, other_cat: other_entry }
                    yield match
                  
def write_comparative_snr(nameA, nameB, tmpA, tmpB):
  sm = SignalMatcher()

  sample_cnt = 0
  for obj in read_object_sequence(tmpA):
    sm.insert(nameA, obj)   

  for obj in read_object_sequence(tmpB):
    sm.insert(nameB, obj)

  snr_sums = {}
  
  all_hours = set()
  all_freqs = set()
  for m in sm.generate_matches():
    sample_cnt += 1
    a = m[nameA]
    b = m[nameB]
    snr_diff = b.snr - a.snr
    unix_ts = truncate_ts_to_hour_unix(a.ts)
    all_hours.add(unix_ts)
    all_freqs.add(a.freq)
    k = (unix_ts, a.freq)
    v = snr_sums.get(k)
    if v is None:
      v = [0, 0.0]
    v[0] += 1
    v[1] += snr_diff

    snr_sums[k] = v


  snr_means = {}
  for k, v in snr_sums.items():
    n = v[0]
    s = v[1]
    result = s/n
    snr_means[k] = result

  del snr_sums

  csv_filename = 'comparative_snr_by_hour.csv'
  all_hours = sorted(list(all_hours))
  all_freqs = sorted(list(all_freqs))
  with open(csv_filename, 'w') as fout:
    fout.write("# input file A: %s\n" % (nameA,))
    fout.write("# input file B: %s\n" % (nameB,))
    fout.write("# hours covered: %d\n" % (len(all_hours),) )
    fout.write("# frequencies covered: %d\n" % (len(all_freqs),) )
    fout.write("# signals received by both antennas: %d\n" % (sample_cnt,))
    fout.write("# ")
    cols = ['hour (unix ts)']
    
    for freq in all_freqs:
      cols.append('%d Hz relative SNR (mean)' % (freq,))
   
    fout.write(','.join(cols))

    fout.write("\n")
    
    for hour in all_hours:
      row = [datetime.datetime.utcfromtimestamp(hour).strftime(DATETIME_FORMAT_FOR_GNUPLOT)]
      for freq in all_freqs:
        k = (hour, freq)
        v = snr_means.get(k)
        if v is None:
          row.append('NaN')
        else:
          row.append('%.1f' % (v,))
      fout.write(','.join(str(x) for x in row))
      fout.write("\n")

  with open('comparative_snr_by_hour.gnuplot', 'w') as fout:
    fout.write("set datafile separator comma\n")
    fout.write("set datafile missing \"NaN\"\n")
    # TODO - put a line at 0 on the Y axis
    filename = "comparative_snr_by_hour.png" 
    fout.write("set terminal pngcairo font \"Arial, 16\" size 1920, 1080\n")
    fout.write("set output \"%s\"\n" % (filename,))
    fout.write("set title \"Average Relative SNR by hour\"\n")
    fout.write("set xlabel \"Hour\"\n")
    fout.write("set ylabel \"SNR above %s\"\n" % (nameA,))

    fout.write("set xdata time\n")
    fout.write(GNUPLOT_TIMEFMT_HRS)
    fout.write(GNUPLOT_X_FORMAT_HRS)
    fout.write("set grid\n")

    accum = []
    for i, freq in enumerate(all_freqs):
      column_idx = 1 + i
      freq_human = "%.3f MHz" % (freq/(10**6.0),)
      # TODO colors that aren't boring
      if i == 0:
        fout.write("plot ")
      else:
        fout.write("\t")

      fout.write("\"%s\" using 1:%d with linespoints lw 2 pt 7 ps 2 lc rgb \"%s\" title \"%s\"" % (csv_filename, column_idx + 1, GNUPLOT_COLORS[i], freq_human,))
      if i != len(all_freqs) - 1:
        fout.write(", \\")
      fout.write("\n")

def filter_signals_distance(rx_gridsquare, station_locator, minimum_distance_km, fin):
  rx_coords = maidenhead.to_location(rx_gridsquare)
  for obj in read_object_sequence(fin):
    callsign = obj.get_tx_station()

    if callsign is None:
      continue
    loc = station_locator.lookup(callsign)

    if loc is None:
      continue # no location known, skip it

    if len(loc) == 4:
      loc = loc + 'mm'
    tx_coords = maidenhead.to_location(loc)
    # TODO switch to pyproj
    dist_km = geopy.distance.great_circle(rx_coords, tx_coords).km

    if dist_km >= minimum_distance_km:
      yield obj

def rx_counts_by_heading(divisions, degrees_per_division, rx_gridsquare, station_locator, minimum_distance_km, fin):
  rx_lat, rx_lon = maidenhead.to_location(rx_gridsquare)
  gds = pyproj.Geod(ellps = 'WGS84')
  
  counts_by_freq = {} 

  for obj in filter_signals_distance(rx_gridsquare, station_locator, minimum_distance_km, fin):
    loc = station_locator.lookup(obj.get_tx_station())
    tx_lat, tx_lon = maidenhead.to_location(loc)
    # TODO - figure out if this is great circle route
    fwd, back, dist = gds.inv(rx_lon, rx_lat, tx_lon, tx_lat)
  
    # interval is [-180.0, 180], add if negative so the range of values is 0-360.0
    if fwd < 0.0:
      fwd += 360.0

    freq = obj.freq 

#    print("%.3f ; %.3f ; %.3f" % (fwd, back, dist,))
    counts = counts_by_freq.get(freq)
    if counts is None:
      counts = { div: 0 for div in range(divisions) } # fill in zeroes
      counts_by_freq[freq] = counts
 
    div = int(math.floor(fwd/degrees_per_division))
    counts[div] += 1

  return counts_by_freq
def write_rx_density_by_heading(rx_gridsquare, station_locator, minimum_distance_km, nameA, nameB, inputA, inputB):
  divisions = 45
  degrees_per_division = 360.0/divisions

  counts_by_freqA = rx_counts_by_heading(divisions, degrees_per_division, rx_gridsquare, station_locator, minimum_distance_km, inputA)
  counts_by_freqB = rx_counts_by_heading(divisions, degrees_per_division, rx_gridsquare, station_locator, minimum_distance_km, inputB)

  x_coords = [i * degrees_per_division for i in range(divisions)]
  x_coords_radians = [ x/180.0 * math.pi for x in x_coords]
  x_coords_radians = x_coords_radians + x_coords_radians
  widths_radians = [degrees_per_division/180.0 * math.pi for _ in range(divisions)]
  widths_radians = widths_radians + widths_radians
  colorsA = [GNUPLOT_COLORS[0] for _ in range(divisions)]
  colorsB = [GNUPLOT_COLORS[1] for _ in range(divisions)]
  all_colors = colorsA + colorsB
  all_freqs = set()
  [all_freqs.add(x) for x in counts_by_freqA.keys()]
  [all_freqs.add(x) for x in counts_by_freqB.keys()]

  for freq in all_freqs: 
    countsA = counts_by_freqA[freq]
    heightsA = [y[1] for y in sorted(countsA.items(), key = lambda x : x[0])]
    countsB = counts_by_freqB[freq]
    heightsB = [y[1] for y in sorted(countsB.items(), key = lambda x : x[0])]
    heights = heightsA + heightsB

    freq_human = '%.3f' % (freq/10**6.0)
    title = "Signals received by heading (%s MHz)" % (freq_human,)
    ax = plt.subplot(111, projection = 'polar', title = title)
    ax.set_title(title)
    
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.bar(x_coords_radians, heights, color = all_colors, width = widths_radians, bottom = 0.0, alpha = 0.44)
    fname = 'rx_count_by_heading_%dHz.png' % (freq,)
    plt.savefig(fname, format = 'png', dpi = 150)

station_locator = StationLocator()
for obj in itertools.chain(*(read_object_sequence(x) for x in (tmpA, tmpB,))):
  callsign = obj.get_tx_station()
  if callsign is not None:    
    station_location = obj.get_maidenhead()
    if station_location is not None:
#      print("%s => %s" % (callsign, station_location,))
      station_locator.insert(callsign, station_location)

rx_gridsquare = 'EM10dk'
minimum_distance_km = 350.0

write_rx_density_by_heading(rx_gridsquare, station_locator, minimum_distance_km, nameA, nameB, tmpA, tmpB)
write_files_for_counts_by_hour(nameA, nameB, tmpA, tmpB)
write_maximum_distance_by_hour(rx_gridsquare, station_locator, nameA, nameB, tmpA, tmpB)
write_comparative_snr(nameA, nameB, tmpA, tmpB)

# Signals RX'd by A, but not by B - by heading - split day/night - exclude groundwave
# Signals RX'd by B, but not by A - by heading - split day/night - exclude groundwave
# probably use this https://matplotlib.org/3.2.2/gallery/pie_and_polar_charts/polar_bar.html#sphx-glr-gallery-pie-and-polar-charts-polar-bar-py

# RX counts / density by heading - split day/night - exclude groundwave
# probably use this https://matplotlib.org/3.2.2/gallery/pie_and_polar_charts/polar_bar.html#sphx-glr-gallery-pie-and-polar-charts-polar-bar-py

# relative SNR by heading & distance - split day/night - exclude groundwave
# maybe we can use this to generate the appropriate plot by using one color for negative and one color positive SNR differences
# https://matplotlib.org/3.2.2/gallery/pie_and_polar_charts/polar_scatter.html#sphx-glr-gallery-pie-and-polar-charts-polar-scatter-py
# the size of the circle indicates the magnitude


tmpA.close()
tmpB.close()
