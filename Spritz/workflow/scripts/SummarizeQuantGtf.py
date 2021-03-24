import sys, os
import numpy as np
import pandas as pd

outtpms = sys.argv[1]
files = sys.argv[2:]
ids = []
tpms = []
for file in files:
    print(f"reading {os.path.basename(file)}")
    currIds = []
    currTpms = []
    for line in open(file):
        if line.startswith("#") : continue
        lins = line.split("\t")
        if lins[2] != "transcript": continue
        info = dict(x.split(' ') for x in lins[8].split('; '))
        transcriptId = info["transcript_id"].strip('"').strip(';')
        covText = info["cov"].strip('"').strip(';')
        fpkmText = info["FPKM"].strip('"').strip(';')
        tpmText = info["TPM"].strip().strip(';').strip('"')
        currIds.append(transcriptId)
        currTpms.append(tpmText)
    if len(ids) == 0: ids = np.asarray(currIds)
    elif all(np.array(currIds) != ids): print("error with ids")
    tpms.append(currTpms)
tpms.insert(0, np.asarray(ids))

print(f"Saving to {outtpms} ...")
dataframe = np.row_stack(tpms)
dataframe[1:,0] = [os.path.basename(file).split(".")[0] for file in files]
pddf = pd.DataFrame(dataframe[1:,1:], index=dataframe[1:,0], columns=dataframe[0,1:]).sort_index()
pddf.T.to_csv(outtpms)
