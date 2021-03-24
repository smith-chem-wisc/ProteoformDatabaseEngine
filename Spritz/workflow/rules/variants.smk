GATK_MEM=24000 # MB
GATK_JAVA=f"--java-options \"-Xmx{GATK_MEM}M -Dsamjdk.compression_level=9\""

rule download_snpeff:
    '''Download and unpack custom SnpEff for annotating variants'''
    output:
        "../resources/SnpEff/snpEff.config",
        "../resources/SnpEff/snpEff.jar",
        filename=temp("../resources/SnpEff_4.3_SmithChemWisc_v2.zip")
    params:
        url="https://github.com/smith-chem-wisc/SnpEff/releases/download/4.3_SCW1/SnpEff_4.3_SmithChemWisc_v2.zip"
    log: "../resources/SnpEffInstall.log"
    conda: "../envs/downloads.yaml"
    shell:
        "(cd ../resources/ && "
        "wget {params.url} && "
        "unzip {output.filename} -d SnpEff) &> {log}"

rule index_fa:
    '''Index genome FASTA file'''
    input: KARYOTYPIC_GENOME_FA
    output: f"{KARYOTYPIC_GENOME_PREFIX}.fa.fai"
    log: f"{KARYOTYPIC_GENOME_PREFIX}.fa.faindex.log"
    conda: "../envs/variants.yaml"
    shell: "samtools faidx {input}"

rule variant_tmpdir:
    output: temp(directory("../resources/tmp")),
    log: "../resources/tmpdir.log"
    conda: "../envs/variants.yaml"
    shell: "mkdir {output} 2> {log}"

rule hisat2_group:
    '''Add read groups to sorted BAM file'''
    input:
        sorted="{dir}/align/combined.sorted.bam",
        tmp="../resources/tmp"
    output:
        grouped=temp("{dir}/variants/combined.sorted.grouped.bam"),
        groupedidx=temp("{dir}/variants/combined.sorted.grouped.bam.bai")
    params:
        gatk_java=GATK_JAVA
    resources:
        mem_mb=GATK_MEM
    log: "{dir}/variants/combined.sorted.grouped.log"
    benchmark: "{dir}/variants/combined.sorted.grouped.benchmark"
    conda: "../envs/variants.yaml"
    shell:
        "(gatk {params.gatk_java} AddOrReplaceReadGroups"
        " -PU platform -PL illumina -SM sample -LB library"
        " -I {input.sorted} -O {output.grouped} --TMP_DIR {input.tmp} && "
        "samtools index {output.grouped}) &> {log}"

rule hisat2_mark:
    '''Mark duplicates in sorted BAM file'''
    input:
        grouped="{dir}/variants/combined.sorted.grouped.bam",
        tmp="../resources/tmp"
    output:
        marked="{dir}/variants/combined.sorted.grouped.marked.bam",
        markedidx="{dir}/variants/combined.sorted.grouped.marked.bam.bai",
        metrics="{dir}/variants/combined.sorted.grouped.marked.metrics"
    params:
        gatk_java=GATK_JAVA
    resources:
        mem_mb=GATK_MEM
    log: "{dir}/variants/combined.sorted.grouped.marked.log"
    benchmark: "{dir}/variants/combined.sorted.grouped.marked.benchmark"
    conda: "../envs/variants.yaml"
    shell:
        "(gatk {params.gatk_java} MarkDuplicates"
        " -I {input.grouped} -O {output.marked} -M {output.metrics}"
        " --ASSUME_SORT_ORDER coordinate --TMP_DIR {input.tmp} && "
        "samtools index {output.marked}) &> {log}"

# Checks if quality encoding is correct, and then splits n cigar reads
rule split_n_cigar_reads:
    '''Check quality scores and split Ns in the cigar reads'''
    input:
        bam="{dir}/variants/combined.sorted.grouped.marked.bam",
        fa=KARYOTYPIC_GENOME_FA,
        fai=f"{KARYOTYPIC_GENOME_PREFIX}.fa.fai",
        fadict=f"{KARYOTYPIC_GENOME_PREFIX}.dict",
        tmp="../resources/tmp"
    output:
        fixed=temp("{dir}/variants/combined.fixedQuals.bam"),
        split=temp("{dir}/variants/combined.sorted.grouped.marked.split.bam"),
        splitidx=temp("{dir}/variants/combined.sorted.grouped.marked.split.bam.bai")
    params:
        gatk_java=GATK_JAVA
    resources:
        mem_mb=GATK_MEM
    log: "{dir}/variants/combined.sorted.grouped.marked.split.log"
    benchmark: "{dir}/variants/combined.sorted.grouped.marked.split.benchmark"
    conda: "../envs/variants.yaml"
    shell:
        "(gatk {params.gatk_java} FixMisencodedBaseQualityReads -I {input.bam} -O {output.fixed} && "
        "gatk {params.gatk_java} SplitNCigarReads -R {input.fa} -I {output.fixed} -O {output.split} --tmp-dir {input.tmp} || " # fix and split
        "gatk {params.gatk_java} SplitNCigarReads -R {input.fa} -I {input.bam} -O {output.split} --tmp-dir {input.tmp}; " # or just split
        "samtools index {output.split}) &> {log}" # always index

rule base_recalibration:
    '''Generate recalibration table and recalibrate BAM file'''
    input:
        knownsites=f"../resources/ensembl/{SPECIES}.ensembl.vcf",
        knownsitesidx=f"../resources/ensembl/{SPECIES}.ensembl.vcf.idx",
        fa=KARYOTYPIC_GENOME_FA,
        bam="{dir}/variants/combined.sorted.grouped.marked.split.bam",
        tmp="../resources/tmp"
    output:
        recaltable=temp("{dir}/variants/combined.sorted.grouped.marked.split.recaltable"),
        recalbam=temp("{dir}/variants/combined.sorted.grouped.marked.split.recal.bam")
    params:
        gatk_java=GATK_JAVA
    resources:
        mem_mb=GATK_MEM
    log: "{dir}/variants/combined.sorted.grouped.marked.split.recal.log"
    benchmark: "{dir}/variants/combined.sorted.grouped.marked.split.recal.benchmark"
    conda: "../envs/variants.yaml"
    shell:
        "(gatk {params.gatk_java} BaseRecalibrator -R {input.fa} -I {input.bam}"
        " --known-sites {input.knownsites} -O {output.recaltable} --tmp-dir {input.tmp} && "
        "gatk {params.gatk_java} ApplyBQSR -R {input.fa} -I {input.bam}"
        " --bqsr-recal-file {output.recaltable} -O {output.recalbam} --tmp-dir {input.tmp} && "
        "samtools index {output.recalbam}) &> {log}"

rule call_gvcf_varaints:
    '''Create genome VCF file'''
    input:
        knownsites=f"../resources/ensembl/{SPECIES}.ensembl.vcf",
        knownsitesidx=f"../resources/ensembl/{SPECIES}.ensembl.vcf.idx",
        fa=KARYOTYPIC_GENOME_FA,
        bam="{dir}/variants/combined.sorted.grouped.marked.split.recal.bam",
        tmp="../resources/tmp"
    output: temp("{dir}/variants/combined.sorted.grouped.marked.split.recal.g.vcf.gz"),
    threads: 8
        # HaplotypeCaller is only fairly efficient with threading;
        # ~14000 regions/min with 24 threads,
        # and ~13000 regions/min with 8 threads,
        # so going with 8 threads max here
    params:
        gatk_java=GATK_JAVA
    resources:
        mem_mb=GATK_MEM
    log: "{dir}/variants/combined.sorted.grouped.marked.split.recal.g.log"
    benchmark: "{dir}/variants/combined.sorted.grouped.marked.split.recal.g.benchmark"
    conda: "../envs/variants.yaml"
    shell:
        "(gatk {params.gatk_java} HaplotypeCaller"
        " --native-pair-hmm-threads {threads}"
        " -R {input.fa} -I {input.bam}"
        " --min-base-quality-score 20 --dont-use-soft-clipped-bases true"
        " --dbsnp {input.knownsites} -O {output} --tmp-dir {input.tmp}"
        " -ERC GVCF --max-mnp-distance 3 && "
        "gatk IndexFeatureFile -I {output}) &> {log}"

rule call_vcf_variants:
    '''Genotype the gVCF for the combined dataset to make VCF'''
    input:
        fa=KARYOTYPIC_GENOME_FA,
        gvcf="{dir}/variants/combined.sorted.grouped.marked.split.recal.g.vcf.gz",
        tmp="../resources/tmp"
    output: "{dir}/variants/combined.sorted.grouped.marked.split.recal.g.gt.vcf" # renamed in next rule
    params:
        gatk_java=GATK_JAVA
    resources:
        mem_mb=GATK_MEM
    log: "{dir}/variants/combined.sorted.grouped.marked.split.recal.g.gt.log"
    benchmark: "{dir}/variants/combined.sorted.grouped.marked.split.recal.g.gt.benchmark"
    conda: "../envs/variants.yaml"
    shell:
        "(gatk {params.gatk_java} GenotypeGVCFs"
        " -R {input.fa} -V {input.gvcf} -O {output} --tmp-dir {input.tmp} && "
        "gatk IndexFeatureFile -I {output}) &> {log}"

rule final_vcf_naming:
    '''Rename VCF to shorter filename'''
    input: "{dir}/variants/combined.sorted.grouped.marked.split.recal.g.gt.vcf"
    output: "{dir}/variants/combined.spritz.vcf"
    log: "{dir}/variants/final_vcf_naming.log"
    conda: "../envs/variants.yaml"
    shell: "mv {input} {output} 2> {log}"

rule variant_annotation_ref:
    '''Generate proteome FASTA and XML for reference database'''
    input:
        f"../resources/SnpEff/data/{REF}/done{REF}.txt",
        snpeff="../resources/SnpEff/snpEff.jar",
        fa=KARYOTYPIC_GENOME_FA,
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
    conda: "../envs/proteogenomics.yaml"
    shell:
        "(java -Xmx{resources.mem_mb}M -jar {input.snpeff} -v -stats {output.html}"
        " -fastaProt {output.protfa} -xmlProt {output.protxml} "
        " {params.ref} {input.vcf}" # no isoforms, with variants
        " > {output.ann}) 2> {log}"

rule variant_annotation_custom:
    input:
        snpeff="../resources/SnpEff/snpEff.jar",
        fa=KARYOTYPIC_GENOME_FA,
        vcf="{dir}/variants/combined.spritz.vcf",
        isoform_reconstruction=[
            "../resources/SnpEff/data/combined.transcripts.genome.gff3/genes.gff",
            "../resources/SnpEff/data/combined.transcripts.genome.gff3/protein.fa",
            "../resources/SnpEff/data/genomes/combined.transcripts.genome.gff3.fa",
            "../resources/SnpEff/data/combined.transcripts.genome.gff3/done.txt"],
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
    conda: "../envs/proteogenomics.yaml"
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
        refprotfa=os.path.join("{dir}/variants/", f"{REF}.{ENSEMBL_VERSION}.protein.fasta"),
        refprotwithdecoysfa=os.path.join("{dir}/variants/", f"{REF}.{ENSEMBL_VERSION}.protein.withdecoys.fasta"),
        refprotwithmodsxml=os.path.join("{dir}/variants/", f"{REF}.{ENSEMBL_VERSION}.protein.withmods.xml.gz"),
    output:
        ann="{dir}/final/combined.spritz.snpeff.vcf",
        protfa="{dir}/final/combined.spritz.snpeff.protein.fasta",
        protwithdecoysfa="{dir}/final/combined.spritz.snpeff.protein.withdecoys.fasta",
        protxmlwithmodsgz="{dir}/final/combined.spritz.snpeff.protein.withmods.xml.gz",
        refprotfa=os.path.join("{dir}/final/", f"{REF}.{ENSEMBL_VERSION}.protein.fasta"),
        refprotwithdecoysfa=os.path.join("{dir}/final/", f"{REF}.{ENSEMBL_VERSION}.protein.withdecoys.fasta"),
        refprotwithmodsxml=os.path.join("{dir}/final/", f"{REF}.{ENSEMBL_VERSION}.protein.withmods.xml.gz"),
    log: "{dir}/variants/finish_isoform_variants.log"
    conda: "../envs/proteogenomics.yaml"
    shell:
        "cp {input.ann} {input.protfa} {input.protwithdecoysfa} {input.protxmlwithmodsgz}"
        " {input.refprotfa} {input.refprotwithdecoysfa} {input.refprotwithmodsxml} {wildcards.dir}/final 2> {log}"

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
    log: "{dir}/variants/finish_isoform_variants.log"
    conda: "../envs/proteogenomics.yaml"
    shell:
        "cp {input.ann} {input.protfa} {input.protwithdecoysfa} {input.protxmlwithmodsgz} {wildcards.dir}/final 2> {log}"
