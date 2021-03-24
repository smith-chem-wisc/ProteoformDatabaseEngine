rule hisat_genome:
    '''Build genome index for hisat2'''
    input:
        fa=KARYOTYPIC_GENOME_FA,
        gtf=ENSEMBL_GFF,
    threads: 12
    output:
        idx=f"{KARYOTYPIC_GENOME_PREFIX}.1.ht2",
        finished=f"../resources/ensembl/done_building_hisat_genome{REF}.txt",
    benchmark: f"../resources/ensembl/{REF}.hisatbuild.benchmark"
    params: prefix=lambda w, input: input.fa[:-3]
    log: f"../resources/ensembl/{REF}.hisatbuild.log"
    conda: "../envs/align.yaml"
    shell:
        "(hisat2-build -p {threads} {input.fa} {params.prefix} && "
        "touch {output.finished}) &> {log}"

rule hisat2_splice_sites:
    '''Fetch the splice sites from the gene model for hisat2'''
    input: ENSEMBL_GFF
    output: f"../resources/ensembl/{SPECIES}.{GENEMODEL_VERSION}.splicesites.txt"
    log: f"../resources/ensembl/{SPECIES}.{GENEMODEL_VERSION}.splicesites.log"
    conda: "../envs/align.yaml"
    shell: "hisat2_extract_splice_sites.py {input} > {output} 2> {log}"

if check('sra'):
    rule prefetch_sras:
        '''Prefetch SRA from GEO SRA'''
        output: temp("{dir}/sra_paired/{sra,[A-Z0-9]+}/{sra}.sra")
        benchmark: "{dir}/{sra}.prefetch.benchmark"
        log: "{dir}/{sra}.log"
        conda: "../envs/sra.yaml"
        shell:
            "prefetch {wildcards.sra}"
            " --output-directory {wildcards.dir}/sra_paired &> {log}"

    rule split_sra: # in the future, could use this to check SE vs PE: https://www.biostars.org/p/139422/
        '''Split fastqs from GEO SRA for quality control and alignment'''
        input: "{dir}/sra_paired/{sra,[A-Z0-9]+}/{sra}.sra"
        output:
            fq1="{dir}/{sra,[A-Z0-9]+}_1.fastq",
            fq2="{dir}/{sra,[A-Z0-9]+}_2.fastq"
        benchmark: "{dir}/{sra}.split.benchmark"
        log: "{dir}/{sra}.log"
        conda: "../envs/sra.yaml"
        shell:
            "fastq-dump -I --outdir {wildcards.dir} --split-files {input} &> {log}"

    rule fastp_sra:
        '''Trim adapters, read quality filtering, make QC outputs'''
        input:
            fq1="{dir}/{sra}_1.fastq",
            fq2="{dir}/{sra}_2.fastq"
        output:
            fq1="{dir}/{sra,[A-Z0-9]+}.trim_1.fastq.gz",
            fq2="{dir}/{sra,[A-Z0-9]+}.trim_2.fastq.gz",
            html="{dir}/{sra,[A-Z0-9]+}.trim.html",
            json="{dir}/{sra,[A-Z0-9]+}.trim.json",
        threads: 6
        log: "{dir}/{sra}.trim.log"
        conda: "../envs/align.yaml"
        params:
            quality=20,
            title="{sra}"
        shell:
            "fastp -q {params.quality}"
            " -i {input.fq1} -I {input.fq2} -o {output.fq1} -O {output.fq2}"
            " -h {output.html} -j {output.json}"
            " -w {threads} -R {params.title} --detect_adapter_for_pe &> {log}"

    rule hisat2_align_bam_sra:
        '''Align trimmed reads'''
        input:
            f"../resources/ensembl/done_building_hisat_genome{REF}.txt",
            fq1="{dir}/{sra}.trim_1.fastq.gz",
            fq2="{dir}/{sra}.trim_2.fastq.gz",
            ss=f"../resources/ensembl/{SPECIES}.{GENEMODEL_VERSION}.splicesites.txt"
        output:
            "{dir}/align/{sra,[A-Z0-9]+}.sra.sorted.bam"
        threads: 12
        params:
            compression="9",
            tempprefix="{dir}/align/{sra}.sra.sorted",
            prefix=KARYOTYPIC_GENOME_PREFIX
        log: "{dir}/align/{sra}.sra.hisat2.log"
        conda: "../envs/align.yaml"
        shell:
            "(hisat2 -p {threads} -x {params.prefix}"
            " -1 {input.fq1} -2 {input.fq2}"
            " --known-splicesite-infile {input.ss} | " # align the suckers
            "samtools view -h -F4 - | " # get mapped reads only
            "samtools sort -l {params.compression} -T {params.tempprefix} -o {output} - && " # sort them
            "samtools index {output}) 2> {log}"

if check('sra_se'):
    rule prefetch_sras_se:
        '''Prefetch SRA from GEO SRA'''
        output: temp("{dir}/sra_single/{sra_se,[A-Z0-9]+}/{sra_se}.sra")
        benchmark: "{dir}/{sra_se}.benchmark"
        log: "{dir}/{sra_se}.log"
        conda: "../envs/sra.yaml"
        shell:
            "prefetch {wildcards.sra_se}"
            " --output-directory {wildcards.dir}/sra_single &> {log}"

    rule split_sras_se:
        input: "{dir}/sra_single/{sra_se,[A-Z0-9]+}/{sra_se}.sra"
        output: "{dir}/{sra_se,[A-Z0-9]+}.fastq" # independent of pe/se
        benchmark: "{dir}/{sra_se}.benchmark"
        log: "{dir}/{sra_se}.log"
        conda: "../envs/sra.yaml"
        shell:
            "fastq-dump -I --outdir {wildcards.dir} --split-files {input} && "
            "mv {wildcards.dir}/{wildcards.sra_se}_1.fastq {output} &> {log}"

    rule fastp_sra_se:
        '''Trim adapters, read quality filtering, make QC outputs'''
        input: "{dir}/{sra_se}.fastq",
        output:
            fq="{dir}/{sra_se,[A-Z0-9]+}.trim.fastq.gz",
            html="{dir}/{sra_se,[A-Z0-9]+}.trim.html",
            json="{dir}/{sra_se,[A-Z0-9]+}.trim.json",
        threads: 6
        log: "{dir}/{sra_se}.trim.log"
        conda: "../envs/align.yaml"
        params:
            quality=20,
            title="{sra_se}"
        shell:
            "fastp -q {params.quality}"
            " -i {input} -o {output.fq}"
            " -h {output.html} -j {output.json}"
            " -w {threads} -R {params.title} --detect_adapter_for_pe &> {log}"

    rule hisat2_align_bam_sra_se:
        '''Align trimmed reads'''
        input:
            f"../resources/ensembl/done_building_hisat_genome{REF}.txt",
            fq="{dir}/{sra_se}.trim.fastq.gz",
            ss=f"../resources/ensembl/{SPECIES}.{GENEMODEL_VERSION}.splicesites.txt"
        output:
            "{dir}/align/{sra_se,[A-Z0-9]+}.sra_se.sorted.bam"
        threads: 12
        params:
            compression="9",
            tempprefix="{dir}/align/{sra_se}.sra_se.sorted",
            prefix=KARYOTYPIC_GENOME_PREFIX
        log: "{dir}/align/{sra_se}.sra_se.hisat2.log"
        conda: "../envs/align.yaml"
        shell:
            "(hisat2 -p {threads} -x {params.prefix}"
            " -U {input.fq}"
            " --known-splicesite-infile {input.ss} | " # align the suckers
            "samtools view -h -F4 - | " # get mapped reads only
            "samtools sort -l {params.compression} -T {params.tempprefix} -o {output} - && " # sort them
            "samtools index {output}) 2> {log}"

if check('fq'):
    rule fastp_fq_uncompressed:
        '''Trim adapters, read quality filtering, make QC outputs'''
        input:
            fq1="{dir}/{fq}_1.fastq",
            fq2="{dir}/{fq}_2.fastq"
        output:
            fq1="{dir}/{fq}.fq.trim_1.fastq.gz",
            fq2="{dir}/{fq}.fq.trim_2.fastq.gz",
            html="{dir}/{fq}.fq.trim.html",
            json="{dir}/{fq}.fq.trim.json",
        threads: 6
        log: "{dir}/{fq}.fq.trim.log"
        conda: "../envs/align.yaml"
        params:
            quality=20,
            title="{fq}"
        shell:
            "fastp -q {params.quality}"
            " -i {input.fq1} -I {input.fq2} -o {output.fq1} -O {output.fq2}"
            " -h {output.html} -j {output.json}"
            " -w {threads} -R {params.title} --detect_adapter_for_pe &> {log}"

    rule fastp_fq:
        '''Trim adapters, read quality filtering, make QC outputs'''
        input:
            fq1="{dir}/{fq}_1.fastq.gz",
            fq2="{dir}/{fq}_2.fastq.gz"
        output:
            fq1="{dir}/{fq}.fq.trim_1.fastq.gz",
            fq2="{dir}/{fq}.fq.trim_2.fastq.gz",
            html="{dir}/{fq}.fq.trim.html",
            json="{dir}/{fq}.fq.trim.json",
        threads: 6
        log: "{dir}/{fq}.fq.trim.log"
        conda: "../envs/align.yaml"
        params:
            quality=20,
            title="{fq}"
        shell:
            "fastp -q {params.quality}"
            " -i {input.fq1} -I {input.fq2} -o {output.fq1} -O {output.fq2}"
            " -h {output.html} -j {output.json}"
            " -w {threads} -R {params.title} --detect_adapter_for_pe &> {log}"

    rule hisat2_align_bam_fq:
        '''Align trimmed reads'''
        input:
            f"../resources/ensembl/done_building_hisat_genome{REF}.txt",
            fq1="{dir}/{fq}.fq.trim_1.fastq.gz",
            fq2="{dir}/{fq}.fq.trim_2.fastq.gz",
            ss=f"../resources/ensembl/{SPECIES}.{GENEMODEL_VERSION}.splicesites.txt"
        output:
            "{dir}/align/{fq}.fq.sorted.bam"
        threads: 12
        params:
            compression="9",
            tempprefix="{dir}/align/{fq}.fq.sorted",
            prefix=KARYOTYPIC_GENOME_PREFIX
        log: "{dir}/align/{fq}.fq.hisat2.log"
        conda: "../envs/align.yaml"
        shell:
            "(hisat2 -p {threads} -x {params.prefix}"
            " -1 {input.fq1} -2 {input.fq2}"
            " --known-splicesite-infile {input.ss} | " # align the suckers
            "samtools view -h -F4 - | " # get mapped reads only
            "samtools sort -l {params.compression} -T {params.tempprefix} -o {output} - && " # sort them
            "samtools index {output}) 2> {log}"

if check('fq_se'):
    rule fastp_fq_se_uncompressed:
        '''Trim adapters, read quality filtering, make QC outputs'''
        input:
            fq1="{dir}/{fq_se}_1.fastq",
        output:
            fq1="{dir}/{fq_se}.fq_se.trim_1.fastq.gz",
            html="{dir}/{fq_se}.fq_se.trim.html",
            json="{dir}/{fq_se}.fq_se.trim.json",
        threads: 6
        log: "{dir}/{fq_se}.fq_se.trim.log"
        conda: "../envs/align.yaml"
        params:
            quality=20,
            title="{fq_se}"
        shell:
            "fastp -q {params.quality}"
            " -i {input.fq1} -o {output.fq1}"
            " -h {output.html} -j {output.json}"
            " -w {threads} -R {params.title} --detect_adapter_for_pe &> {log}"

    rule fastp_fq_se:
        '''Trim adapters, read quality filtering, make QC outputs'''
        input:
            fq1="{dir}/{fq_se}_1.fastq.gz",
        output:
            fq1="{dir}/{fq_se}.fq_se.trim_1.fastq.gz",
            html="{dir}/{fq_se}.fq_se.trim.html",
            json="{dir}/{fq_se}.fq_se.trim.json",
        threads: 6
        log: "{dir}/{fq_se}.fq_se.trim.log"
        conda: "../envs/align.yaml"
        params:
            quality=20,
            title="{fq_se}"
        shell:
            "fastp -q {params.quality}"
            " -i {input.fq1} -o {output.fq1}"
            " -h {output.html} -j {output.json}"
            " -w {threads} -R {params.title} --detect_adapter_for_pe &> {log}"

    rule hisat2_align_bam_fq_se:
        '''Align trimmed reads'''
        input:
            f"../resources/ensembl/done_building_hisat_genome{REF}.txt",
            fq1="{dir}/{fq_se}.fq_se.trim_1.fastq.gz",
            ss=f"../resources/ensembl/{SPECIES}.{GENEMODEL_VERSION}.splicesites.txt"
        output:
            "{dir}/align/{fq_se}.fq_se.sorted.bam"
        threads: 12
        params:
            compression="9",
            tempprefix="{dir}/align/{fq_se}.fq_se.sorted",
            prefix=KARYOTYPIC_GENOME_PREFIX
        log: "{dir}/align/{fq_se}.fq_se.hisat2.log"
        conda: "../envs/align.yaml"
        shell:
            "(hisat2 -p {threads} -x {params.prefix}"
            " -U {input.fq1}"
            " --known-splicesite-infile {input.ss} | " # align the suckers
            "samtools view -h -F4 - | " # get mapped reads only
            "samtools sort -l {params.compression} -T {params.tempprefix} -o {output} - && " # sort them
            "samtools index {output}) 2> {log}"

rule hisat2_merge_bams:
    '''Merge the BAM files for each sample'''
    input:
        lambda w:
            ([] if not check('sra') else expand("{{dir}}/align/{sra}.sra.sorted.bam", sra=config["sra"])) + \
            ([] if not check('sra_se') else expand("{{dir}}/align/{sra_se}.sra_se.sorted.bam", sra_se=config["sra_se"])) + \
            ([] if not check('fq') else expand("{{dir}}/align/{fq}.fq.sorted.bam", fq=config["fq"])) + \
            ([] if not check('fq_se') else expand("{{dir}}/align/{fq_se}.fq_se.sorted.bam", fq_se=config["fq_se"])),
    output:
        sorted="{dir}/align/combined.sorted.bam",
        stats="{dir}/align/combined.sorted.stats"
    params:
        compression="9",
        tempprefix=lambda w, input: os.path.splitext(input[0])[0]
    log: "{dir}/align/combined.sorted.log"
    conda: "../envs/align.yaml"
    threads: 12
    resources: mem_mb=16000
    shell:
        "(ls {input} | "
        "{{ read firstbam; "
        "samtools view -h ""$firstbam""; "
        "while read bam; do samtools view ""$bam""; done; }} | "
        "samtools view -ubS - | "
        "samtools sort -@ {threads} -l {params.compression} -T {params.tempprefix} -o {output.sorted} - && "
        "samtools index {output.sorted} && "
        "samtools flagstat -@ {threads} {output.sorted} > {output.stats}) 2> {log}"
