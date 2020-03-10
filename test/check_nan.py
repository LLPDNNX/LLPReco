import sys
import os

import numpy as np
import ROOT
import pandas as pd
import uproot
from array import array
import csv

branches = ["global_*", "csv_*", "cpf_*", "npf_*", "sv_*", "length_*", "muon_*", "electron_*", "jetorigin_*"]

file = "nano.root"
features = uproot.open(file)["Events"].arrays(branches).keys()

for i in range(5):
    try:
        tree = uproot.open(file)["Events"]
        break
    except:
        time.sleep(10)

nan = 0
inf = 0

for feature in features:

    #print(feature)
    feature_array = tree.array(feature).flatten()

    if np.isnan(feature_array).any():
        print("NAN in " + str(feature))
        nan += 1
    else:
        pass

    if np.isinf(feature_array).any():
        print("INF in " + str(feature))
        inf += 1
    else:
        pass

if (nan + inf > 0):
    print("number of NAN + INF: " + str(nan + inf))
    sys.exit(1)
