SPECIES_LOWER = config['species'].lower()
PROTOCOL = "http"

rule download_ensembl_references:
    output:
        gfa=f"data/ensembl/{REF}.dna.primary_assembly.fa",
        gff3=f"data/ensembl/{REF}.{config['release']}.gff3",
        pfa=f"data/ensembl/{REF}.pep.all.fa",
    params:
        primary=f"{PROTOCOL}://ftp.ensembl.org/pub/release-{ENSEMBL_VERSION}//fasta/{SPECIES_LOWER}/dna/{REF}.dna.primary_assembly.fa.gz",
        toplevel=f"{PROTOCOL}://ftp.ensembl.org/pub/release-{ENSEMBL_VERSION}//fasta/{SPECIES_LOWER}/dna/{REF}.dna.toplevel.fa.gz",
        gff=f"{PROTOCOL}://ftp.ensembl.org/pub/release-{ENSEMBL_VERSION}/gff3/{SPECIES_LOWER}/{REF}.{ENSEMBL_VERSION}.gff3.gz",
        pep=f"{PROTOCOL}://ftp.ensembl.org/pub/release-{ENSEMBL_VERSION}//fasta/{SPECIES_LOWER}/pep/{REF}.pep.all.fa.gz",
    benchmark: "data/ensembl/downloads.benchmark"
    log: "data/ensembl/downloads.log"
    conda: "environments/downloads.yaml"
    shell:
        "((wget -O - {params.primary} || wget -O - {params.toplevel}) | gunzip -c - > {output.gfa} && "
        "wget -O - {params.gff} | gunzip -c - > {output.gff3} && "
        "wget -O - {params.pep} | gunzip -c - > {output.pfa}) 2> {log}"

if SPECIES_LOWER == "homo_sapiens":
    rule download_dbsnp_vcf:
        '''Download dbsnp known variant sites if we are analyzing human data'''
        input: f"ChromosomeMappings/{config['genome']}_UCSC2ensembl.txt"
        output: f"data/ensembl/{config['species']}.ensembl.vcf",
        params:
            vcf="https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/common_all_20180418.vcf.gz"
        benchmark: "data/ensembl/downloads_dbsnp_vcf.benchmark"
        log: "data/ensembl/downloads_dbsnp_vcf.log"
        conda: "environments/downloads.yaml"
        shell:
            "(wget -O - {params.vcf} | zcat - | python scripts/convert_ucsc2ensembl.py > {output}) 2> {log}"
else:
    rule download_ensembl_vcf:
        '''
        Use Ensembl known variant sites if we are analyzing nonhuman data.
        Note that Ensembl has started listing variants for each chromosome for human, but not other species, but that may change
        '''
        output: f"data/ensembl/{config['species']}.ensembl.vcf",
        params:
            vcf1 = f"http://ftp.ensembl.org/pub/release-{ENSEMBL_VERSION}/variation/vcf/{SPECIES_LOWER}/{SPECIES}.vcf.gz",
            vcf2 = f"http://ftp.ensembl.org/pub/release-{ENSEMBL_VERSION}/variation/vcf/{SPECIES_LOWER}/{SPECIES_LOWER}.vcf.gz",
        benchmark: "data/ensembl/downloads_ensembl_vcf.benchmark"
        log: "data/ensembl/downloads_ensembl_vcf.log"
        conda: "environments/downloads.yaml"
        shell: "((wget -O - {params.vcf1} || wget -O - {params.vcf2}) | zcat - | python scripts/clean_vcf.py > {output}) 2> {log}"

rule index_ensembl_vcf:
    input: f"data/ensembl/{config['species']}.ensembl.vcf"
    output: f"data/ensembl/{config['species']}.ensembl.vcf.idx"
    log: f"data/ensembl/{config['species']}.ensembl.vcf.idx.log"
    conda: "environments/variants.yaml"
    shell: "gatk IndexFeatureFile -I {input} 2> {log}"

rule download_chromosome_mappings:
    output: f"ChromosomeMappings/{config['genome']}_UCSC2ensembl.txt"
    params: url="https://github.com/dpryan79/ChromosomeMappings.git"
    log: "data/download_chromosome_mappings.log"
    conda: "environments/downloads.yaml"
    shell:
        "(if [ -d ChromosomeMappings ]; then rm -rf ChromosomeMappings; fi && "
        "git clone {params.url}) 2> {log}"

rule reorder_genome_fasta:
    input: f"data/ensembl/{REF}.dna.primary_assembly.fa"
    output: f"data/ensembl/{REF}.dna.primary_assembly.karyotypic.fa"
    benchmark: "data/ensembl/karyotypic_order.benchmark"
    log: "data/ensembl/karyotypic_order.log"
    conda: "environments/downloads.yaml"
    shell: "python scripts/karyotypic_order.py 2> {log}"

rule dict_fa:
    input: f"data/ensembl/{config['species']}.{config['genome']}.dna.primary_assembly.karyotypic.fa"
    output: f"data/ensembl/{config['species']}.{config['genome']}.dna.primary_assembly.karyotypic.dict"
    log: f"data/ensembl/{config['species']}.{config['genome']}.dna.primary_assembly.karyotypic.dict.log"
    conda: "environments/variants.yaml"
    shell: "gatk CreateSequenceDictionary -R {input} -O {output} 2> {log}"

rule tmpdir:
    output:
        temp(directory("tmp")),
        temp(directory("temporary")),
    log: "data/tmpdir.log"
    shell:
        "mkdir tmp temporary 2> {log}"
