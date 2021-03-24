import os

# Variables used by many of the rules
MIN_SPRITZ_VERSION = "0.3.0" # should be the same here, common.smk, and MainWindow.xml.cs
SPECIES = config["species"]
GENOME_VERSION = config["genome"]
ENSEMBL_VERSION = config["release"]
GENEMODEL_VERSION = f"{GENOME_VERSION}.{ENSEMBL_VERSION}"
REF = f"{SPECIES}.{GENOME_VERSION}"
GENOME_FA = f"../resources/ensembl/{REF}.dna.primary_assembly.fa"
KARYOTYPIC_GENOME_PREFIX = f"../resources/ensembl/{REF}.dna.primary_assembly.karyotypic"
KARYOTYPIC_GENOME_FA = f"{KARYOTYPIC_GENOME_PREFIX}.fa"
ENSEMBL_GFF = f"../resources/ensembl/{SPECIES}.{GENEMODEL_VERSION}.gff3"
TEST_GENOME_FA = f"../resources/ensembl/202122.fa"
TEST_ENSEMBL_GFF = f"../resources/ensembl/202122.gff3"
FA=GENOME_FA # for analysis; can also be TEST_GENOME_FA
GFF3=ENSEMBL_GFF # for analysis; can also be TEST_ENSEMBL_GFF
REFSTAR_PREFIX = f"../resources/ensembl/{SPECIES}.{GENEMODEL_VERSION}RsemStar/RsemStarReference"
REFSTAR_FOLDER = f"../resources/ensembl/{SPECIES}.{GENEMODEL_VERSION}RsemStar/"
REF_PREFIX = f"../resources/ensembl/{SPECIES}.{GENEMODEL_VERSION}Rsem/RsemReference"
REF_FOLDER = f"../resources/ensembl/{SPECIES}.{GENEMODEL_VERSION}Rsem/"
UNIPROTXML=f"../resources/uniprot/{config['species']}.protein.xml.gz" #"../resources/Homo_sapiens_202022.xml.gz"
UNIPROTFASTA=f"../resources/uniprot/{config['species']}.protein.fasta" #"../resources/Homo_sapiens_202022.xml.gz"
TRANSFER_MOD_DLL="../TransferUniProtModifications/TransferUniProtModifications/bin/Release/netcoreapp3.1/TransferUniProtModifications.dll"

def all_output(wildcards):
    '''Gets the final output files depending on the configuration'''
    outputs = ["prose.txt"]
    min_version = MIN_SPRITZ_VERSION.split('.')
    this_version = config['spritz_version'].split('.')
    if not "spritz_version" in config or any([min_version[i] > this_version[i] for i in range(len(min_version))]):
        outputs = ["please_update_spritz.txt"]
    else:
        outputs.append(os.path.join("variants/", f"done{REF}.{ENSEMBL_VERSION}.txt")) # reference
        if "variant" in config["analyses"]:
            outputs.append("final/combined.spritz.snpeff.protein.withmods.xml.gz") # variants
        if "isoform" in config["analyses"]:
            outputs.append("final/combined.spritz.isoform.protein.withmods.xml.gz") # isoforms
        if "quant" in config["analyses"]:
            outputs.extend([
                "final/transcript_reference_quant.tpms.csv",
                "final/gene_reference_quant.tpms.csv"]) # reference quant with stringtie
        if "variant" in config["analyses"] and "isoform" in config["analyses"]:
            outputs.append("final/combined.spritz.isoformvariants.protein.withmods.xml.gz") # isoform variants
        if "isoform" in config["analyses"] and "quant" in config["analyses"]:
            outputs.extend([
                "final/transcript_custom_quant.tpms.csv",
                "final/gene_custom_quant.tpms.csv"]) # isoform quant with stringtie
    expanded_outputs = expand(
        [os.path.join("{dir}", file) for file in outputs],
        dir=config["analysis_directory"])
    return expanded_outputs

def check(field):
    '''Checks whether or not a field is contained in the configuration'''
    return field in config and config[field] is not None and len(config[field]) > 0
