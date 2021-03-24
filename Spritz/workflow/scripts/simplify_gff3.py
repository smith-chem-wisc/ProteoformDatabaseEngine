import sys
gff = sys.argv[1]
for line in open(gff):
    linesplit = line.split('\t')
    emptyOrComment = len(linesplit) < 3 or line.startswith("#")
    if not emptyOrComment and not linesplit[2] == "exon" and not linesplit[2].endswith("UTR"): sys.stdout.write(line)
