GATK_MEM=16000 # MB
GATK_JAVA=f"--java-options \"-Xmx{GATK_MEM}M -Dsamjdk.compression_level=9\""
REF=config["species"] + "." + config["genome"]
GENOME_VERSION = config["genome"]

rule download_snpeff:
    output: "SnpEff/snpEff.config", "SnpEff/snpEff.jar", temp("SnpEff_4.3_SmithChemWisc_v2.zip")
    log: "data/SnpEffInstall.log"
    shell:
        "(wget https://github.com/smith-chem-wisc/SnpEff/releases/download/4.3_SCW1/SnpEff_4.3_SmithChemWisc_v2.zip && "
        "unzip SnpEff_4.3_SmithChemWisc_v2.zip -d SnpEff) &> {log}"

rule index_fa:
    input: "data/ensembl/" + REF + ".dna.primary_assembly.karyotypic.fa"
    output: "data/ensembl/" + REF + ".dna.primary_assembly.karyotypic.fa.fai"
    shell: "samtools faidx {input}"

rule hisat2_groupmark_bam:
    input:
        sorted="{dir}/align/combined.sorted.bam",
        tmp=directory("tmp")
    output:
        grouped=temp("{dir}/variants/combined.sorted.grouped.bam"),
        groupedidx=temp("{dir}/variants/combined.sorted.grouped.bam.bai"),
        marked="{dir}/variants/combined.sorted.grouped.marked.bam",
        markedidx="{dir}/variants/combined.sorted.grouped.marked.bam.bai",
        metrics="{dir}/variants/combined.sorted.grouped.marked.metrics"
    resources: mem_mb=GATK_MEM
    log: "{dir}/variants/combined.sorted.grouped.marked.log"
    benchmark: "{dir}/variants/combined.sorted.grouped.marked.benchmark"
    shell:
        "(gatk {GATK_JAVA} AddOrReplaceReadGroups -PU platform  -PL illumina -SM sample -LB library -I {input.sorted} -O {output.grouped} -SO coordinate --TMP_DIR {input.tmp} && "
        "samtools index {output.grouped} && "
        "gatk {GATK_JAVA} MarkDuplicates -I {output.grouped} -O {output.marked} -M {output.metrics} --TMP_DIR {input.tmp} -AS true && "
        "samtools index {output.marked}) &> {log}"

# Checks if quality encoding is correct, and then splits n cigar reads
rule split_n_cigar_reads:
    input:
        bam="{dir}/variants/combined.sorted.grouped.marked.bam",
        fa="data/ensembl/" + REF + ".dna.primary_assembly.karyotypic.fa",
        fai="data/ensembl/" + REF + ".dna.primary_assembly.karyotypic.fa.fai",
        fadict="data/ensembl/" + REF + ".dna.primary_assembly.karyotypic.dict",
        tmp=directory("tmp")
    output:
        fixed=temp("{dir}/variants/combined.fixedQuals.bam"),
        split=temp("{dir}/variants/combined.sorted.grouped.marked.split.bam"),
        splitidx=temp("{dir}/variants/combined.sorted.grouped.marked.split.bam.bai")
    resources: mem_mb=GATK_MEM
    log: "{dir}/variants/combined.sorted.grouped.marked.split.log"
    benchmark: "{dir}/variants/combined.sorted.grouped.marked.split.benchmark"
    shell:
        "(gatk {GATK_JAVA} FixMisencodedBaseQualityReads -I {input.bam} -O {output.fixed} && "
        "gatk {GATK_JAVA} SplitNCigarReads -R {input.fa} -I {output.fixed} -O {output.split} --tmp-dir {input.tmp} || " # fix and split
        "gatk {GATK_JAVA} SplitNCigarReads -R {input.fa} -I {input.bam} -O {output.split} --tmp-dir {input.tmp}; " # or just split
        "samtools index {output.split}) &> {log}" # always index

rule base_recalibration:
    input:
        knownsites="data/ensembl/" + config["species"] + ".ensembl.vcf",
        knownsitesidx="data/ensembl/" + config["species"] + ".ensembl.vcf.idx",
        fa="data/ensembl/" + REF + ".dna.primary_assembly.karyotypic.fa",
        bam="{dir}/variants/combined.sorted.grouped.marked.split.bam",
        tmp=directory("tmp")
    output:
        recaltable=temp("{dir}/variants/combined.sorted.grouped.marked.split.recaltable"),
        recalbam=temp("{dir}/variants/combined.sorted.grouped.marked.split.recal.bam")
    resources: mem_mb=GATK_MEM
    log: "{dir}/variants/combined.sorted.grouped.marked.split.recal.log"
    benchmark: "{dir}/variants/combined.sorted.grouped.marked.split.recal.benchmark"
    shell:
        "(gatk {GATK_JAVA} BaseRecalibrator -R {input.fa} -I {input.bam} --known-sites {input.knownsites} -O {output.recaltable} --tmp-dir {input.tmp} && "
        "gatk {GATK_JAVA} ApplyBQSR -R {input.fa} -I {input.bam} --bqsr-recal-file {output.recaltable} -O {output.recalbam} --tmp-dir {input.tmp} && "
        "samtools index {output.recalbam}) &> {log}"

rule call_gvcf_varaints:
    input:
        knownsites="data/ensembl/" + config["species"] + ".ensembl.vcf",
        knownsitesidx="data/ensembl/" + config["species"] + ".ensembl.vcf.idx",
        fa="data/ensembl/" + REF + ".dna.primary_assembly.karyotypic.fa",
        bam="{dir}/variants/combined.sorted.grouped.marked.split.recal.bam",
        tmp=directory("tmp")
    output: temp("{dir}/variants/combined.sorted.grouped.marked.split.recal.g.vcf.gz"),
    threads: 8
        # HaplotypeCaller is only fairly efficient with threading;
        # ~14000 regions/min with 24 threads,
        # and ~13000 regions/min with 8 threads,
        # so going with 8 threads max here
    resources: mem_mb=GATK_MEM
    log: "{dir}/variants/combined.sorted.grouped.marked.split.recal.g.log"
    benchmark: "{dir}/variants/combined.sorted.grouped.marked.split.recal.g.benchmark"
    shell:
        "(gatk {GATK_JAVA} HaplotypeCaller"
        " --native-pair-hmm-threads {threads}"
        " -R {input.fa} -I {input.bam}"
        " --min-base-quality-score 20 --dont-use-soft-clipped-bases true"
        " --dbsnp {input.knownsites} -O {output} --tmp-dir {input.tmp}"
        " -ERC GVCF --max-mnp-distance 3 && "
        "gatk IndexFeatureFile -F {output}) &> {log}"

rule call_vcf_variants:
    input:
        fa="data/ensembl/" + REF + ".dna.primary_assembly.karyotypic.fa",
        gvcf="{dir}/variants/combined.sorted.grouped.marked.split.recal.g.vcf.gz",
        tmp=directory("tmp")
    output: "{dir}/variants/combined.sorted.grouped.marked.split.recal.g.gt.vcf" # renamed in next rule
    resources: mem_mb=GATK_MEM
    log: "{dir}/variants/combined.sorted.grouped.marked.split.recal.g.gt.log"
    benchmark: "{dir}/variants/combined.sorted.grouped.marked.split.recal.g.gt.benchmark"
    shell:
        "(gatk {GATK_JAVA} GenotypeGVCFs -R {input.fa} -V {input.gvcf} -O {output} --tmp-dir {input.tmp} && "
        "gatk IndexFeatureFile -F {output}) &> {log}"

rule final_vcf_naming:
    input: "{dir}/variants/combined.sorted.grouped.marked.split.recal.g.gt.vcf"
    output: "{dir}/variants/combined.spritz.vcf"
    shell: "mv {input} {output}"

rule filter_indels:
    input:
        fa="data/ensembl/" + REF + ".dna.primary_assembly.karyotypic.fa",
        vcf="{dir}/variants/combined.spritz.vcf"
    output: "{dir}/variants/combined.spritz.noindels.vcf"
    log: "{dir}/variants/combined.spritz.noindels.log"
    benchmark: "{dir}/variants/combined.spritz.noindels.benchmark"
    shell:
        "(gatk SelectVariants --select-type-to-exclude INDEL -R {input.fa} -V {input.vcf} -O {output} && "
        "gatk IndexFeatureFile -F {output}) &> {log}"

rule variant_annotation_ref:
    input:
        "SnpEff/data/" + REF + "/done" + REF + ".txt",
        snpeff="SnpEff/snpEff.jar",
        fa="data/ensembl/" + REF + ".dna.primary_assembly.karyotypic.fa",
        vcf="{dir}/variants/combined.spritz.vcf",
    output:
        ann="{dir}/variants/combined.spritz.snpeff.vcf",
        html="{dir}/variants/combined.spritz.snpeff.html",
        genesummary="{dir}/variants/combined.spritz.snpeff.genes.txt",
        protfa="{dir}/variants/combined.spritz.snpeff.protein.fasta",
        protxml="{dir}/variants/combined.spritz.snpeff.protein.xml"
    params: ref=REF, # no isoform reconstruction
    resources: mem_mb=16000
    log: "{dir}/variants/combined.spritz.snpeff.log"
    benchmark: "{dir}/variants/combined.spritz.snpeff.benchmark"
    shell:
        "(java -Xmx{resources.mem_mb}M -jar {input.snpeff} -v -stats {output.html}"
        " -fastaProt {output.protfa} -xmlProt {output.protxml} "
        " {params.ref} {input.vcf}" # no isoforms, with variants
        " > {output.ann}) 2> {log}"

rule variant_annotation_custom:
    input:
        snpeff="SnpEff/snpEff.jar",
        fa="data/ensembl/" + REF + ".dna.primary_assembly.karyotypic.fa",
        vcf="{dir}/variants/combined.spritz.vcf",
        isoform_reconstruction=[
            "SnpEff/data/combined.transcripts.genome.gff3/genes.gff",
            "SnpEff/data/combined.transcripts.genome.gff3/protein.fa",
            "SnpEff/data/genomes/combined.transcripts.genome.gff3.fa",
            "SnpEff/data/combined.transcripts.genome.gff3/done.txt"],
    output:
        ann="{dir}/variants/combined.spritz.isoformvariants.vcf",
        html="{dir}/variants/combined.spritz.isoformvariants.html",
        genesummary="{dir}/variants/combined.spritz.isoformvariants.genes.txt",
        protfa="{dir}/variants/combined.spritz.isoformvariants.protein.fasta",
        protxml=temp("{dir}/variants/combined.spritz.isoformvariants.protein.xml"),
    params: ref="combined.transcripts.genome.gff3" # with isoforms
    resources: mem_mb=GATK_MEM
    log: "{dir}/variants/combined.spritz.isoformvariants.log"
    benchmark: "{dir}/variants/combined.spritz.isoformvariants.benchmark"
    shell:
        "(java -Xmx{resources.mem_mb}M -jar {input.snpeff} -v -stats {output.html}"
        " -fastaProt {output.protfa} -xmlProt {output.protxml}"
        " {params.ref} {input.vcf}" # with isoforms and variants
        " > {output.ann}) 2> {log}"

rule finish_variants:
    '''Copy final output files from variant workflow to main directory'''
    input:
        ann="{dir}/variants/combined.spritz.snpeff.vcf",
        protfa="{dir}/variants/combined.spritz.snpeff.protein.fasta",
        protwithdecoysfa="{dir}/variants/combined.spritz.snpeff.protein.withdecoys.fasta",
        protxmlwithmodsgz="{dir}/variants/combined.spritz.snpeff.protein.withmods.xml.gz",
        refprotfa="{dir}/variants/" + REF + "." + ENSEMBL_VERSION + ".protein.fasta",
        refprotwithdecoysfa="{dir}/variants/" + REF + "." + ENSEMBL_VERSION + ".protein.withdecoys.fasta",
        refprotwithmodsxml="{dir}/variants/" + REF + "." + ENSEMBL_VERSION + ".protein.withmods.xml.gz",
    output:
        ann="{dir}/final/combined.spritz.snpeff.vcf",
        protfa="{dir}/final/combined.spritz.snpeff.protein.fasta",
        protwithdecoysfa="{dir}/final/combined.spritz.snpeff.protein.withdecoys.fasta",
        protxmlwithmodsgz="{dir}/final/combined.spritz.snpeff.protein.withmods.xml.gz",
        refprotfa="{dir}/final/" + REF + "." + ENSEMBL_VERSION + ".protein.fasta",
        refprotwithdecoysfa="{dir}/final/" + REF + "." + ENSEMBL_VERSION + ".protein.withdecoys.fasta",
        refprotwithmodsxml="{dir}/final/" + REF + "." + ENSEMBL_VERSION + ".protein.withmods.xml.gz",
    shell:
        "cp {input.ann} {input.protfa} {input.protwithdecoysfa} {input.protxmlwithmodsgz}"
        " {input.refprotfa} {input.refprotwithdecoysfa} {input.refprotwithmodsxml} {wildcards.dir}/final"

rule finish_isoform_variants:
    '''Copy final output files from isoform-variant workflow to main directory'''
    input:
        ann="{dir}/variants/combined.spritz.isoformvariants.vcf",
        protfa="{dir}/variants/combined.spritz.isoformvariants.protein.fasta",
        protwithdecoysfa="{dir}/variants/combined.spritz.isoformvariants.protein.withdecoys.fasta",
        protxmlwithmodsgz="{dir}/variants/combined.spritz.isoformvariants.protein.withmods.xml.gz",
    output:
        ann="{dir}/final/combined.spritz.isoformvariants.vcf",
        protfa="{dir}/final/combined.spritz.isoformvariants.protein.fasta",
        protwithdecoysfa="{dir}/final/combined.spritz.isoformvariants.protein.withdecoys.fasta",
        protxmlwithmodsgz="{dir}/final/combined.spritz.isoformvariants.protein.withmods.xml.gz",
    shell:
        "cp {input.ann} {input.protfa} {input.protwithdecoysfa} {input.protxmlwithmodsgz} {wildcards.dir}/final"
