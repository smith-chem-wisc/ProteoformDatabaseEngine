from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import yaml

with open("config/config.yaml", 'r') as stream:
   data = yaml.safe_load(stream)

species = data["species"]
version = data["genome"]

fasta_sequences = SeqIO.parse(open(f"../resources/ensembl/{species}.{version}.dna.primary_assembly.fa"), 'fasta')

with open("../resources/ensembl/202122.fa", "w") as out:
    for seq_record in fasta_sequences:
        if seq_record.id.startswith("20") or seq_record.id.startswith("21") or seq_record.id.startswith("22"):
            out.write(seq_record.format("fasta"))
