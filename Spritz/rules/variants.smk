GATK_MEM=16000 # MB
GATK_JAVA=f"--java-options \"-Xmx{GATK_MEM}M -Dsamjdk.compression_level=9\""

rule download_snpeff:
    output:
        "SnpEff/snpEff.config",
        "SnpEff/snpEff.jar",
        filename=temp("SnpEff_4.3_SmithChemWisc_v2.zip")
    params:
        url="https://github.com/smith-chem-wisc/SnpEff/releases/download/4.3_SCW1/SnpEff_4.3_SmithChemWisc_v2.zip"
    log: "data/SnpEffInstall.log"
    conda: "environments/downloads.yaml"
    shell:
        "(wget {params.url} && unzip {output.filename} -d SnpEff) &> {log}"

rule index_fa:
    input: f"data/ensembl/{REF}.dna.primary_assembly.karyotypic.fa"
    output: f"data/ensembl/{REF}.dna.primary_assembly.karyotypic.fa.fai"
    log: f"data/ensembl/{REF}.dna.primary_assembly.karyotypic.fa.log"
    conda: "environments/variants.yaml"
    shell: "samtools faidx {input}"

rule hisat2_groupmark_bam:
    input:
        sorted="{dir}/align/combined.sorted.bam",
        tmp="tmp"
    output:
        grouped=temp("{dir}/variants/combined.sorted.grouped.bam"),
        groupedidx=temp("{dir}/variants/combined.sorted.grouped.bam.bai"),
        marked="{dir}/variants/combined.sorted.grouped.marked.bam",
        markedidx="{dir}/variants/combined.sorted.grouped.marked.bam.bai",
        metrics="{dir}/variants/combined.sorted.grouped.marked.metrics"
    params:
        gatk_java=GATK_JAVA
    resources:
        mem_mb=GATK_MEM
    log: "{dir}/variants/combined.sorted.grouped.marked.log"
    benchmark: "{dir}/variants/combined.sorted.grouped.marked.benchmark"
    conda: "environments/variants.yaml"
    shell:
        "(gatk {params.gatk_java} AddOrReplaceReadGroups -PU platform  -PL illumina -SM sample -LB library -I {input.sorted} -O {output.grouped} -SO coordinate --TMP_DIR {input.tmp} && "
        "samtools index {output.grouped} && "
        "gatk {params.gatk_java} MarkDuplicates -I {output.grouped} -O {output.marked} -M {output.metrics} --TMP_DIR {input.tmp} -AS true && "
        "samtools index {output.marked}) &> {log}"

# Checks if quality encoding is correct, and then splits n cigar reads
rule split_n_cigar_reads:
    input:
        bam="{dir}/variants/combined.sorted.grouped.marked.bam",
        fa=f"data/ensembl/{REF}.dna.primary_assembly.karyotypic.fa",
        fai=f"data/ensembl/{REF}.dna.primary_assembly.karyotypic.fa.fai",
        fadict=f"data/ensembl/{REF}.dna.primary_assembly.karyotypic.dict",
        tmp="tmp"
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
    conda: "environments/variants.yaml"
    shell:
        "(gatk {params.gatk_java} FixMisencodedBaseQualityReads -I {input.bam} -O {output.fixed} && "
        "gatk {params.gatk_java} SplitNCigarReads -R {input.fa} -I {output.fixed} -O {output.split} --tmp-dir {input.tmp} || " # fix and split
        "gatk {params.gatk_java} SplitNCigarReads -R {input.fa} -I {input.bam} -O {output.split} --tmp-dir {input.tmp}; " # or just split
        "samtools index {output.split}) &> {log}" # always index

rule base_recalibration:
    input:
        knownsites=f"data/ensembl/{config['species']}.ensembl.vcf",
        knownsitesidx=f"data/ensembl/{config['species']}.ensembl.vcf.idx",
        fa=f"data/ensembl/{REF}.dna.primary_assembly.karyotypic.fa",
        bam="{dir}/variants/combined.sorted.grouped.marked.split.bam",
        tmp="tmp"
    output:
        recaltable=temp("{dir}/variants/combined.sorted.grouped.marked.split.recaltable"),
        recalbam=temp("{dir}/variants/combined.sorted.grouped.marked.split.recal.bam")
    params:
        gatk_java=GATK_JAVA
    resources:
        mem_mb=GATK_MEM
    log: "{dir}/variants/combined.sorted.grouped.marked.split.recal.log"
    benchmark: "{dir}/variants/combined.sorted.grouped.marked.split.recal.benchmark"
    conda: "environments/variants.yaml"
    shell:
        "(gatk {params.gatk_java} BaseRecalibrator -R {input.fa} -I {input.bam} --known-sites {input.knownsites} -O {output.recaltable} --tmp-dir {input.tmp} && "
        "gatk {params.gatk_java} ApplyBQSR -R {input.fa} -I {input.bam} --bqsr-recal-file {output.recaltable} -O {output.recalbam} --tmp-dir {input.tmp} && "
        "samtools index {output.recalbam}) &> {log}"

rule call_gvcf_varaints:
    input:
        knownsites=f"data/ensembl/{config['species']}.ensembl.vcf",
        knownsitesidx=f"data/ensembl/{config['species']}.ensembl.vcf.idx",
        fa=f"data/ensembl/{REF}.dna.primary_assembly.karyotypic.fa",
        bam="{dir}/variants/combined.sorted.grouped.marked.split.recal.bam",
        tmp="tmp"
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
    conda: "environments/variants.yaml"
    shell:
        "(gatk {params.gatk_java} HaplotypeCaller"
        " --native-pair-hmm-threads {threads}"
        " -R {input.fa} -I {input.bam}"
        " --min-base-quality-score 20 --dont-use-soft-clipped-bases true"
        " --dbsnp {input.knownsites} -O {output} --tmp-dir {input.tmp}"
        " -ERC GVCF --max-mnp-distance 3 && "
        "gatk IndexFeatureFile -I {output}) &> {log}"

rule call_vcf_variants:
    input:
        fa=f"data/ensembl/{REF}.dna.primary_assembly.karyotypic.fa",
        gvcf="{dir}/variants/combined.sorted.grouped.marked.split.recal.g.vcf.gz",
        tmp="tmp"
    output: "{dir}/variants/combined.sorted.grouped.marked.split.recal.g.gt.vcf" # renamed in next rule
    params:
        gatk_java=GATK_JAVA
    resources:
        mem_mb=GATK_MEM
    log: "{dir}/variants/combined.sorted.grouped.marked.split.recal.g.gt.log"
    benchmark: "{dir}/variants/combined.sorted.grouped.marked.split.recal.g.gt.benchmark"
    conda: "environments/variants.yaml"
    shell:
        "(gatk {params.gatk_java} GenotypeGVCFs -R {input.fa} -V {input.gvcf} -O {output} --tmp-dir {input.tmp} && "
        "gatk IndexFeatureFile -I {output}) &> {log}"

rule final_vcf_naming:
    input: "{dir}/variants/combined.sorted.grouped.marked.split.recal.g.gt.vcf"
    output: "{dir}/variants/combined.spritz.vcf"
    log: "{dir}/variants/final_vcf_naming.log"
    shell: "mv {input} {output} 2> {log}"

rule filter_indels:
    input:
        fa=f"data/ensembl/{REF}.dna.primary_assembly.karyotypic.fa",
        vcf="{dir}/variants/combined.spritz.vcf"
    output: "{dir}/variants/combined.spritz.noindels.vcf"
    log: "{dir}/variants/combined.spritz.noindels.log"
    benchmark: "{dir}/variants/combined.spritz.noindels.benchmark"
    conda: "environments/variants.yaml"
    shell:
        "(gatk SelectVariants --select-type-to-exclude INDEL -R {input.fa} -V {input.vcf} -O {output} && "
        "gatk IndexFeatureFile -I {output}) &> {log}"

rule variant_annotation_ref:
    input:
        f"SnpEff/data/{REF}/done{REF}.txt",
        snpeff="SnpEff/snpEff.jar",
        fa=f"data/ensembl/{REF}.dna.primary_assembly.karyotypic.fa",
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
    conda: "environments/variants.yaml"
    shell:
        "(java -Xmx{resources.mem_mb}M -jar {input.snpeff} -v -stats {output.html}"
        " -fastaProt {output.protfa} -xmlProt {output.protxml} "
        " {params.ref} {input.vcf}" # no isoforms, with variants
        " > {output.ann}) 2> {log}"

rule variant_annotation_custom:
    input:
        snpeff="SnpEff/snpEff.jar",
        fa=f"data/ensembl/{REF}.dna.primary_assembly.karyotypic.fa",
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
    conda: "environments/variants.yaml"
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
    shell:
        "cp {input.ann} {input.protfa} {input.protwithdecoysfa} {input.protxmlwithmodsgz} {wildcards.dir}/final 2> {log}"
