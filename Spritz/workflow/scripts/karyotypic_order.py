from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import yaml

with open("config/config.yaml", 'r') as stream:
   data = yaml.safe_load(stream)

species = data["species"]
version = data["genome"]

fa = f"../resources/ensembl/{species}.{version}.dna.primary_assembly.fa"

ordered = []
chrn = [str(x) for x in range(1, 23)]
chrs = []
x = ""
y = ""
m = ""
gl = []
ki = []
other = []
for seq in SeqIO.parse(open(fa), 'fasta'):
    if seq.id.split(" ")[0] in chrn: chrs.append(seq)
    elif seq.id.split(" ")[0].startswith("X"): x = seq
    elif seq.id.split(" ")[0].startswith("Y"): y = seq
    elif seq.id.split(" ")[0].startswith("MT"): m = seq
    elif seq.id.split(" ")[0].startswith("GL"): gl.append(seq)
    elif seq.id.split(" ")[0].startswith("KI"): ki.append(seq)
    else: other.append(seq)
chrs.sort(key = lambda seq: int(seq.id.split(" ")[0]))
ordered.extend(chrs)
ordered.append(x)
ordered.append(y)
ordered.append(m)
ordered.extend(gl)
ordered.extend(ki)
ordered.extend(other)

with open(f"../resources/ensembl/{species}.{version}.dna.primary_assembly.karyotypic.fa", "w") as out:
    for seq in ordered:
        out.write(seq.format("fasta"))
