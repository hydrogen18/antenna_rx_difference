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

from astropy.coordinates import EarthLocation
import astropy.units as u
from astroplan import Observer

class SignalObservation(object):
  __slots__ = ('ts', 'freq', 'direction', 'mode', 'rssi', 'freq_offset', 'msg',)
  def __init__(self, ts, freq, direction, mode, rssi, _, freq_offset, msg):
    self.ts = ts
    self.freq = freq
    self.direction = direction
    self.mode = mode
    self.rssi = rssi
    self.freq_offset = freq_offset
    self.msg = msg

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

    ts, freq, direction, mode, rssi, _, freq_offset, msg = m.groups()

    dot = freq.find('.')
    freq_scale = 1
    if dot != -1:
      dot_dist = len(freq) - dot
      freq_scale = 10**(dot_dist-1)
     

    freq = int(freq.replace('.', '')) * freq_scale
    freq_offset = int(freq_offset)
    rssi = int(rssi)
    
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
    yield SignalObservation(ts, freq, direction, mode, rssi, None, freq_offset, msg)


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
    truncated_ts = datetime.datetime(obj.ts.year, obj.ts.month, obj.ts.day, obj.ts.hour)
    unix_ts = int(time.mktime(truncated_ts.timetuple()))
    k = (unix_ts, obj.freq,)
    v = data.get(k)
    if v is None:
      v = 0
    data[k] = v + 1
  return data

def write_files_for_counts_by_hour(nameA, nameB, tmpA, tmpB):
  dataA = count_by_freq_and_hour(tmpA)
  dataB = count_by_freq_and_hour(tmpB)
  all_freqs = set()
  [all_freqs.add(x[1]) for x in dataA.keys()]
  [all_freqs.add(x[1]) for x in dataB.keys()]

  all_hours = set()
  [all_hours.add(x[0]) for x in dataA.keys()]
  [all_hours.add(x[0]) for x in dataB.keys()]
  all_hours = sorted(list(all_hours))

  freq_ordering = sorted(list(all_freqs))
  csv_filename = 'counts_by_hour.csv'
  with open(csv_filename, 'w') as fout:
    fout.write("# input file A: %s\n" % (nameA,))
    fout.write("# input file B: %s\n" % (nameB,))
    fout.write("# hours covered: %d\n" % (len(all_hours),) )
    fout.write("# frequencies covered: %d\n" % (len(all_freqs),) )
    fout.write("# ")
    cols = ['hour (unix ts)']
    
    for freq in freq_ordering:
      cols.append('%d Hz count (%s)' % (freq, nameA,))
      cols.append('%d Hz count (%s)' % (freq, nameB,))
   
    fout.write(','.join(cols))

    fout.write("\n")

    for hour in all_hours:
      row = [datetime.datetime.utcfromtimestamp(hour).strftime('%Y/%m/%d %H')]
      for freq in freq_ordering:
        k = (hour, freq,)
        v_a = dataA.get(k, 0)
        v_b = dataB.get(k, 0)
        row.append(v_a)
        row.append(v_b)

      fout.write(','.join(str(x) for x in row))
      fout.write("\n")

  with open('counts_by_hour.gnuplot', 'w') as fout:
    fout.write("set datafile separator comma\n")
    for i, freq in enumerate(freq_ordering):
      filename = "counts_by_hour_%dHz.png" % (freq,)
      freq_human = "%.3f MHz" % (freq/(10**6.0),)

      fout.write("set terminal pngcairo font \"Arial, 16\" size 1920, 1080\n")
      fout.write("set output \"%s\"\n" % (filename,))

      fout.write("set title \"Signals received per hour (%s)\"\n" % (freq_human,))
      fout.write("set xlabel \"Hour\"\n")
      fout.write("set ylabel \"Count\"\n")
      # TODO legend
      # TODO x axis as hours
      a_column_idx = 1 + (i*2)
      b_column_idx = 1 + (i*2 + 1)

      fout.write("set xdata time\n")
      fout.write("set timefmt \"%Y/%m/%d %H\"\n")
#      fout.write("set xrange [\"2020/08/01 00\":\"2020/08/31 23\"]\n")
      fout.write("set format x \"%m-%d\\n%H\"\n")
      fout.write("set grid\n")

      # TODO colors
      fout.write("plot \"%s\" using 1:%d with linespoints lw 3 pt 7 ps 2 lc rgb \"black\" title \"%s\", \\\n" % (csv_filename, a_column_idx + 1, nameA,))
      fout.write("     \"%s\" using 1:%d with linespoints lw 3 pt 7 ps 2 lc rgb \"red\" title \"%s\"\n" % (csv_filename, b_column_idx + 1, nameB,))
      fout.write("\n")

    
write_files_for_counts_by_hour(nameA, nameB, tmpA, tmpB)

tmpA.close()
tmpB.close()
