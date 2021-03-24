import sys, yaml

with open("config/config.yaml", 'r') as stream:
   data = yaml.safe_load(stream)

species = data["species"]
version = data["genome"]

ucsc2ensembl={}
for line in open(f"../resources/ChromosomeMappings/{version}_UCSC2ensembl.txt"):
    linesplit=line.strip().split("\t")
    if len(linesplit) <= 1: continue
    ucsc2ensembl[linesplit[0]] = linesplit[1]

with open(f"../resources/ensembl/{species}.ensembl.vcf","w") as ensembl:
    chrs={}
    max_chr=0
    for line in sys.stdin:
        # header
        if line.startswith("#"):
            sys.stdout.write(line)
            continue

        # change chr from UCSC to Ensembl
        splitline = line.split("\t")
        if len(splitline) > 1 and splitline[0] in ucsc2ensembl:
            splitline[0] = ucsc2ensembl[splitline[0]]
        if splitline[0].isdigit() and int(splitline[0]) > max_chr:
            max_chr = int(splitline[0])
        lineline="\t".join(splitline)
        if not splitline[0] in chrs:
            chrs[splitline[0]] = [lineline]
        else:
            chrs[splitline[0]].append(lineline)

    # order and output
    ordered = []
    chrn = [str(x) for x in range(1, max_chr + 1)]
    chrn.extend(["X","Y","MT"])
    chrn.extend([ccc for ccc in chrs.keys() if ccc.startswith("GL")])
    chrn.extend([ccc for ccc in chrs.keys() if ccc.startswith("KI")])
    otherchr = [ccc for ccc in chrs.keys() if ccc not in chrn]
    chrn.extend(otherchr)
    for chr in chrn:
        if chr in chrs:
            for line in chrs[chr]:
                sys.stdout.write(line)
