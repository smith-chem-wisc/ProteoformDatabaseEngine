import re, sys

bases = set(['A','C','T','G'])
for line in sys.stdin:
    # header
    if line.startswith("#"):
        sys.stdout.write(line)
        continue

    # remove any lines with empty alleles
    splitline = line.split("\t")
    if '' in splitline[0:7]: continue

    # remove unparsable alleles
    if not any(set(splitline[4]).intersection(bases)): continue

    # remove duplicate alleles
    duplicate = re.findall("([a-zA-z]\.)", splitline[4])
    if len(duplicate) > 1: continue

    sys.stdout.write(line)
