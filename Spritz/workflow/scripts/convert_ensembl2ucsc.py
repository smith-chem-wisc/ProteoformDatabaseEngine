import sys, os

file, newfile = sys.argv[1:]

ensembl2ucsc={}
for line in open("../resources/ChromosomeMappings/GRCh38_ensembl2UCSC.txt"):
    linesplit=line.strip().split("\t")
    if len(linesplit) <= 1: continue
    ensembl2ucsc[linesplit[0]] = linesplit[1]

with open(newfile,"w") as ucsc:
    chrs={}
    max_chr=0
    for line in open(file):
        # header
        if line.startswith("#"):
            ucsc.write(line)
            continue

        # change chr from UCSC to Ensembl
        splitline = line.split("\t")
        if '' in splitline[0:7]: continue
        if splitline[0].isdigit() and int(splitline[0]) > max_chr:
            max_chr = int(splitline[0])
        if len(splitline) > 1 and splitline[0] in ensembl2ucsc:
            splitline[0] = ensembl2ucsc[splitline[0]]
        lineline="\t".join(splitline)
        if not splitline[0] in chrs:
            chrs[splitline[0]] = [lineline]
        else:
            chrs[splitline[0]].append(lineline)

    # order and output
    ordered = []
    chrn = [f"chr{x}" for x in range(1, max_chr + 1)]
    chrn.extend(["chrX", "chrY", "chrM"])
    # chrn.extend([ccc for ccc in chrs.keys() if ccc.startswith("GL")])
    # chrn.extend([ccc for ccc in chrs.keys() if ccc.startswith("KI")])
    otherchr = [ccc for ccc in chrs.keys() if ccc not in chrn]
    chrn.extend(otherchr)
    for chr in chrn:
        if chr in chrs:
            for line in chrs[chr]:
                ucsc.write(line)
