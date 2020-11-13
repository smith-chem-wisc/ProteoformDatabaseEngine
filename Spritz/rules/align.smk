rule hisat_genome:
    '''Build genome index for hisat2'''
    input:
        fa=f"data/ensembl/{REF}.dna.primary_assembly.karyotypic.fa",
        gtf=f"data/ensembl/{REF}.{config['release']}.gff3",
    threads: 12
    output:
        idx=f"data/ensembl/{REF}.dna.primary_assembly.karyotypic.1.ht2",
        finished=f"data/ensembl/done_building_hisat_genome{REF}.txt",
    benchmark: f"data/ensembl/{REF}.hisatbuild.benchmark"
    params: ref=REF
    log: f"data/ensembl/{REF}.hisatbuild.log"
    shell:
        "(hisat2-build -p {threads} data/ensembl/{params.ref}.dna.primary_assembly.karyotypic.fa"
        " data/ensembl/{params.ref}.dna.primary_assembly.karyotypic && touch {output.finished}) &> {log}"

rule hisat2_splice_sites:
    '''Fetch the splice sites from the gene model for hisat2'''
    input: f"data/ensembl/{REF}.{config['release']}.gff3"
    output: f"data/ensembl/{REF}.{config['release']}.splicesites.txt"
    log: f"data/ensembl/{REF}.{config['release']}.splicesites.log"
    shell: "hisat2_extract_splice_sites.py {input} > {output} 2> {log}"

if check('sra'):
    rule download_sras: # in the future, could use this to check SE vs PE: https://www.biostars.org/p/139422/
        '''Download fastqs from GEO SRA for quality control and alignment'''
        output:
            fq1="{dir}/{sra,[A-Z0-9]+}_1.fastq",
            fq2="{dir}/{sra,[A-Z0-9]+}_2.fastq"
        benchmark: "{dir}/{sra}.benchmark"
        log: "{dir}/{sra}.log"
        threads: 4
        shell:
            "fasterq-dump -b 10MB -c 100MB -m 1000MB -p --threads {threads}" # use 10x the default memory allocation for larger SRAs
            " --split-files --temp {wildcards.dir} --outdir {wildcards.dir} {wildcards.sra} 2> {log}"

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
        params:
            quality=20,
            title="{sra}"
        shell:
            "fastp -q {params.quality} "
            "-i {input.fq1} -I {input.fq2} -o {output.fq1} -O {output.fq2} "
            "-h {output.html} -j {output.json} "
            "-w {threads} -R {params.title} --detect_adapter_for_pe &> {log}"

    rule hisat2_align_bam_sra:
        '''Align trimmed reads'''
        input:
            f"data/ensembl/{REF}.dna.primary_assembly.karyotypic.1.ht2",
            fq1="{dir}/{sra}.trim_1.fastq.gz",
            fq2="{dir}/{sra}.trim_2.fastq.gz",
            ss=f"data/ensembl/{REF}.{config['release']}.splicesites.txt"
        output:
            "{dir}/align/{sra,[A-Z0-9]+}.sra.sorted.bam"
        threads: 12
        params:
            compression="9",
            tempprefix="{dir}/align/{sra}.sra.sorted",
            ref=REF
        log: "{dir}/align/{sra}.sra.hisat2.log"
        shell:
            "(hisat2 -p {threads} -x data/ensembl/{params.ref}.dna.primary_assembly.karyotypic "
            "-1 {input.fq1} -2 {input.fq2} "
            "--known-splicesite-infile {input.ss} | " # align the suckers
            "samtools view -h -F4 - | " # get mapped reads only
            "samtools sort -l {params.compression} -T {params.tempprefix} -o {output} -) 2> {log} && " # sort them
            "samtools index {output}"

if check('sra_se'):
    rule download_sras_se:
        output: "{dir}/{sra_se,[A-Z0-9]+}.fastq" # independent of pe/se
        benchmark: "{dir}/{sra_se}.benchmark"
        log: "{dir}/{sra_se}.log"
        threads: 4
        shell:
            "fasterq-dump -b 10MB -c 100MB -m 1000MB -p --threads {threads}" # use 10x the default memory allocation for larger SRAs
            " --split-files --temp {wildcards.dir} --outdir {wildcards.dir} {wildcards.sra_se} 2> {log}"

    rule fastp_sra_se:
        '''Trim adapters, read quality filtering, make QC outputs'''
        input: "{dir}/{sra_se}.fastq",
        output:
            fq="{dir}/{sra_se,[A-Z0-9]+}.trim.fastq.gz",
            html="{dir}/{sra_se,[A-Z0-9]+}.trim.html",
            json="{dir}/{sra_se,[A-Z0-9]+}.trim.json",
        threads: 6
        log: "{dir}/{sra_se}.trim.log"
        params:
            quality=20,
            title="{sra_se}"
        shell:
            "fastp -q {params.quality} "
            "-i {input} -o {output.fq} "
            "-h {output.html} -j {output.json} "
            "-w {threads} -R {params.title} --detect_adapter_for_pe &> {log}"

    rule hisat2_align_bam_sra_se:
        '''Align trimmed reads'''
        input:
            f"data/ensembl/{REF}.dna.primary_assembly.karyotypic.1.ht2",
            fq="{dir}/{sra_se}.trim.fastq.gz",
            ss=f"data/ensembl/{REF}.{config['release']}.splicesites.txt"
        output:
            "{dir}/align/{sra_se,[A-Z0-9]+}.sra_se.sorted.bam"
        threads: 12
        params:
            compression="9",
            tempprefix="{dir}/align/{sra_se}.sra_se.sorted",
            ref=REF
        log: "{dir}/align/{sra_se}.sra_se.hisat2.log"
        shell:
            "(hisat2 -p {threads} -x data/ensembl/{params.ref}.dna.primary_assembly.karyotypic "
            "-U {input.fq} "
            "--known-splicesite-infile {input.ss} | " # align the suckers
            "samtools view -h -F4 - | " # get mapped reads only
            "samtools sort -l {params.compression} -T {params.tempprefix} -o {output} -) 2> {log} && " # sort them
            "samtools index {output}"

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
        params:
            quality=20,
            title="{fq}"
        shell:
            "fastp -q {params.quality} "
            "-i {input.fq1} -I {input.fq2} -o {output.fq1} -O {output.fq2} "
            "-h {output.html} -j {output.json} "
            "-w {threads} -R {params.title} --detect_adapter_for_pe &> {log}"

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
        params:
            quality=20,
            title="{fq}"
        shell:
            "fastp -q {params.quality} "
            "-i {input.fq1} -I {input.fq2} -o {output.fq1} -O {output.fq2} "
            "-h {output.html} -j {output.json} "
            "-w {threads} -R {params.title} --detect_adapter_for_pe &> {log}"

    rule hisat2_align_bam_fq:
        '''Align trimmed reads'''
        input:
            f"data/ensembl/{REF}.dna.primary_assembly.karyotypic.1.ht2",
            fq1="{dir}/{fq}.fq.trim_1.fastq.gz",
            fq2="{dir}/{fq}.fq.trim_2.fastq.gz",
            ss=f"data/ensembl/{REF}.{config['release']}.splicesites.txt"
        output:
            "{dir}/align/{fq}.fq.sorted.bam"
        threads: 12
        params:
            compression="9",
            tempprefix="{dir}/align/{fq}.fq.sorted",
            ref=REF
        log: "{dir}/align/{fq}.fq.hisat2.log"
        shell:
            "(hisat2 -p {threads} -x data/ensembl/{params.ref}.dna.primary_assembly.karyotypic "
            "-1 {input.fq1} -2 {input.fq2} "
            "--known-splicesite-infile {input.ss} | " # align the suckers
            "samtools view -h -F4 - | " # get mapped reads only
            "samtools sort -l {params.compression} -T {params.tempprefix} -o {output} -) 2> {log} && " # sort them
            "samtools index {output}"

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
        params:
            quality=20,
            title="{fq_se}"
        shell:
            "fastp -q {params.quality} "
            "-i {input.fq1} -o {output.fq1} "
            "-h {output.html} -j {output.json} "
            "-w {threads} -R {params.title} --detect_adapter_for_pe &> {log}"

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
        params:
            quality=20,
            title="{fq_se}"
        shell:
            "fastp -q {params.quality} "
            "-i {input.fq1} -o {output.fq1} "
            "-h {output.html} -j {output.json} "
            "-w {threads} -R {params.title} --detect_adapter_for_pe &> {log}"

    rule hisat2_align_bam_fq_se:
        '''Align trimmed reads'''
        input:
            f"data/ensembl/{REF}.dna.primary_assembly.karyotypic.1.ht2",
            fq1="{dir}/{fq_se}.fq_se.trim_1.fastq.gz",
            ss=f"data/ensembl/{REF}.{config['release']}.splicesites.txt"
        output:
            "{dir}/align/{fq_se}.fq_se.sorted.bam"
        threads: 12
        params:
            compression="9",
            tempprefix="{dir}/align/{fq_se}.fq_se.sorted",
            ref=REF
        log: "{dir}/align/{fq_se}.fq_se.hisat2.log"
        shell:
            "(hisat2 -p {threads} -x data/ensembl/{params.ref}.dna.primary_assembly.karyotypic "
            "-U {input.fq1} "
            "--known-splicesite-infile {input.ss} | " # align the suckers
            "samtools view -h -F4 - | " # get mapped reads only
            "samtools sort -l {params.compression} -T {params.tempprefix} -o {output} -) 2> {log} && " # sort them
            "samtools index {output}"

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
