from joblib import Parallel, delayed
import multiprocessing

from rnnf import Rnn
from qrnnf import Rnn as Rnnq
import h5py
import argparse
import os
import datetime
import numpy as np
from extract_events import extract_events

def scale94(X):
  m25 = np.percentile(X[:,0], 25)
  m75 = np.percentile(X[:,0], 75)
  s50 = np.median(X[:,2])
  me25 = -0.3
  me75= 0.3
  se50 = 0.6103758
  ret = np.array(X)
  scale = (me75 - me25) / (m75 - m25)
  m25 *= scale
  shift = me25 - m25
  ret[:,0] = X[:,0] * scale + shift
  ret[:,1] = ret[:,0]**2
  
  sscale = se50 / s50

  ret[:,2] = X[:,2] * sscale
  return ret

def scale(X):
  m25 = np.percentile(X[:,0], 25)
  m75 = np.percentile(X[:,0], 75)
  s50 = np.median(X[:,2])
  me25 = 0.07499809
  me75 = 0.26622871
  se50 = 0.6103758
  ret = np.array(X)
  scale = (me75 - me25) / (m75 - m25)
  m25 *= scale
  shift = me25 - m25
  ret[:,0] = X[:,0] * scale + shift
  ret[:,1] = ret[:,0]**2
  
  sscale = se50 / s50

  ret[:,2] = X[:,2] * sscale
  return ret

def get_events(h5):
  if not args.event_detect:
    try:
      e = h5["Analyses/Basecall_RNN_1D_000/BaseCalled_template/Events"]
      return e
    except:
      pass
    try:
      e = h5["Analyses/Basecall_1D_000/BaseCalled_template/Events"]
      return e
    except:
      pass

  return extract_events(h5, args.chemistry)

def basecall(i, filename):
  o_file_name = "res_R9/output" + str(i) + ".fasta"
  
  if os.path.isfile(o_file_name) and os.path.getsize(o_file_name) > 0: # already basecalled
    return 0

  try:
    h5 = h5py.File(filename, "r")
    events = get_events(h5)
    if events is None:
      print "No events in file %s" % filename
      h5.close()
      return 0

    if len(events) < 300:
      print "Read %s too short, not basecalling" % filename
      h5.close()
      return 0

    events = events[50:-50][:args.max_events]
    mean = events["mean"]
    std = events["stdv"]
    length = events["length"]
    scale_f = scale if args.chemistry == 'r9' else scale94
    X = np.array(np.vstack([mean, mean*mean, std, length]).T, dtype=np.float32)

    if len(X) < 2500 or args.cut_data == False:
      X = scale_f(X)
      o1, o2 = ntwk.predict(X)
    else:
      preds1 = []
      preds2 = []
      for i in range(0, len(X), 2000):
        o1, o2 = ntwk.predict(scale_f(X[i:i+2500]))
        if i > 0:
          o1 = o1[250:]
          o2 = o2[250:]
        if i + 2500 < len(X):
          o1 = o1[:-250]
          o2 = o2[:-250]
        preds1.append(o1)
        preds2.append(o2)

      o1 = np.vstack(preds1)
      o2 = np.vstack(preds2)

    o1m = (np.argmax(o1, 1))
    o2m = (np.argmax(o2, 1))

    om = np.vstack((o1m,o2m)).reshape((-1,),order='F')
    output = "".join(map(lambda x: alph[x], om)).replace("N", "")
    
    output_file = open("res_R9/output" + str(i) + ".fasta", "w")
    print >>output_file, ">%s_template_deepnano" % filename
    print >>output_file, output
    output_file.flush()
    output_file.close()
    h5.close()
    return len(events)
  except Exception as e:
    print "Read %s failed with %s" % (filename, e)
    return 0

alph = "ACGTN"

parser = argparse.ArgumentParser()
parser.add_argument('--chemistry', choices=['r9', 'r9.4'], default='r9.4')
parser.add_argument('--output', type=str, default="output.fasta")
parser.add_argument('--directory', type=str, default='', help="Directory where read files are stored")
parser.add_argument('--watch', type=str, default='', help='Watched directory')
parser.add_argument('reads', type=str, nargs='*')
parser.add_argument('--debug', dest='debug', action='store_true')
parser.add_argument('--no-debug', dest='debug', action='store_false')
parser.add_argument('--event-detect', dest='event_detect', action='store_true')
parser.add_argument('--max-events', type=int, default=50000, help='Max amount of events to basecall')
parser.add_argument('--cut-data', dest='cut_data', action='store_true', help='Cut data into smaller chunks and basecall them separatelly. Helps in case of bad preprocessing.')
parser.set_defaults(debug=False)
parser.set_defaults(event_detect=False)
parser.set_defaults(cut_data=False)

args = parser.parse_args()

assert len(args.reads) != 0 or len(args.directory) != 0 or len(args.watch) != 0, "Nothing to basecall"

ntwks = {"r9": os.path.join("networks", "r9.pkl"), "r9.4": os.path.join("networks", "r94.pkl")}

ntwk = Rnn() if args.chemistry == 'r9' else Rnnq()
ntwk.load(ntwks[args.chemistry])


if len(args.reads) or len(args.directory) != 0:
  files = args.reads
  if len(args.directory):
    files += [os.path.join(args.directory, x) for x in os.listdir(args.directory)]  
  Parallel(n_jobs=6)(delayed(basecall)(i, read) for i, read in enumerate(files))
