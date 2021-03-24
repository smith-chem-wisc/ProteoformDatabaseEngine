if check('sra'):
    rule quantify_transcripts_sra:
        input:
            bam="{dir}/align/{sra}.sra.sorted.bam",
            gff="{dir}/isoforms/combined.gtf",
        output:
            ["{dir}/isoforms/{sra}_quant/t_data.ctab",
            "{dir}/isoforms/{sra}_quant/e_data.ctab",
            "{dir}/isoforms/{sra}_quant/i_data.ctab",
            "{dir}/isoforms/{sra}_quant/e2t.ctab",
            "{dir}/isoforms/{sra}_quant/i2t.ctab"],
            genes="{dir}/isoforms/{sra}_quant/{sra}.sra.gene.quant.tab",
            gtf=temp("{dir}/isoforms/{sra}_quant/{sra}.sra.transcript.quant.gtf"),
            gtfgz="{dir}/isoforms/{sra}_quant/{sra}.sra.transcript.quant.gtf.gz"
        threads: 4
        log: "{dir}/isoforms/{sra}_quant/{sra}.sra.quant.gtf.log"
        params: outfolder=lambda w, output: os.path.dirname(output.genes)
        benchmark: "{dir}/isoforms/{sra}_quant/{sra}.sra.quant.gtf.benchmark"
        conda: "../envs/quant.yaml"
        shell:
            "stringtie {input.bam} -p {threads} -G {input.gff} -A {output.genes} -o {output.gtf} -eb {params.outfolder} && "
            "gzip -k {output.gtf} 2> {log}"

    rule quantify_ref_transcripts_sra:
        input:
            bam="{dir}/align/{sra}.sra.sorted.bam",
            gff=GFF3,
        output:
            ["{dir}/isoforms/{sra}_quant_ref/t_data.ctab",
            "{dir}/isoforms/{sra}_quant_ref/e_data.ctab",
            "{dir}/isoforms/{sra}_quant_ref/i_data.ctab",
            "{dir}/isoforms/{sra}_quant_ref/e2t.ctab",
            "{dir}/isoforms/{sra}_quant_ref/i2t.ctab"],
            genes="{dir}/isoforms/{sra}_quant_ref/{sra}.sra.gene.quant_ref.tab",
            gtf=temp("{dir}/isoforms/{sra}_quant_ref/{sra}.sra.transcript.quant_ref.gtf"),
            gtfgz="{dir}/isoforms/{sra}_quant_ref/{sra}.sra.transcript.quant_ref.gtf.gz"
        threads: 4
        log: "{dir}/isoforms/{sra}_quant_ref/{sra}.sra.quant_ref.gtf.log"
        params: outfolder=lambda w, output: os.path.dirname(output.genes)
        benchmark: "{dir}/isoforms/{sra}_quant_ref/{sra}.sra.quant_ref.gtf.benchmark"
        conda: "../envs/quant.yaml"
        shell:
            "stringtie {input.bam} -p {threads} -G {input.gff} -A {output.genes} -o {output.gtf} -eb {params.outfolder} && "
            "gzip -k {output.gtf} 2> {log}"

if check('sra_se'):
    rule quantify_transcripts_sra_se:
        input:
            bam="{dir}/align/{sra_se}.sra_se.sorted.bam",
            gff="{dir}/isoforms/combined.gtf",
        output:
            ["{dir}/isoforms/{sra_se}_quant/t_data.ctab",
            "{dir}/isoforms/{sra_se}_quant/e_data.ctab",
            "{dir}/isoforms/{sra_se}_quant/i_data.ctab",
            "{dir}/isoforms/{sra_se}_quant/e2t.ctab",
            "{dir}/isoforms/{sra_se}_quant/i2t.ctab"],
            genes="{dir}/isoforms/{sra_se}_quant/{sra_se}.sra_se.gene.quant.tab",
            gtf=temp("{dir}/isoforms/{sra_se}_quant/{sra_se}.sra_se.transcript.quant.gtf"),
            gtfgz="{dir}/isoforms/{sra_se}_quant/{sra_se}.sra_se.transcript.quant.gtf.gz"
        threads: 4
        log: "{dir}/isoforms/{sra_se}_quant/{sra_se}.sra_se.quant.gtf.log"
        params: outfolder=lambda w, output: os.path.dirname(output.genes)
        benchmark: "{dir}/isoforms/{sra_se}_quant/{sra_se}.sra_se.quant.gtf.benchmark"
        conda: "../envs/quant.yaml"
        shell:
            "stringtie {input.bam} -p {threads} -G {input.gff} -A {output.genes} -o {output.gtf} -eb {params.outfolder} && "
            "gzip -k {output.gtf} 2> {log}"

    rule quantify_ref_transcripts_sra_se:
        input:
            bam="{dir}/align/{sra_se}.sra_se.sorted.bam",
            gff=GFF3,
        output:
            ["{dir}/isoforms/{sra_se}_quant_ref/t_data.ctab",
            "{dir}/isoforms/{sra_se}_quant_ref/e_data.ctab",
            "{dir}/isoforms/{sra_se}_quant_ref/i_data.ctab",
            "{dir}/isoforms/{sra_se}_quant_ref/e2t.ctab",
            "{dir}/isoforms/{sra_se}_quant_ref/i2t.ctab"],
            genes="{dir}/isoforms/{sra_se}_quant_ref/{sra_se}.sra_se.gene.quant_ref.tab",
            gtf=temp("{dir}/isoforms/{sra_se}_quant_ref/{sra_se}.sra_se.transcript.quant_ref.gtf"),
            gtfgz="{dir}/isoforms/{sra_se}_quant_ref/{sra_se}.sra_se.transcript.quant_ref.gtf.gz"
        threads: 4
        log: "{dir}/isoforms/{sra_se}_quant_ref/{sra_se}.sra_se.quant_ref.gtf.log"
        params: outfolder=lambda w, output: os.path.dirname(output.genes)
        benchmark: "{dir}/isoforms/{sra_se}_quant_ref/{sra_se}.sra_se.quant_ref.gtf.benchmark"
        conda: "../envs/quant.yaml"
        shell:
            "stringtie {input.bam} -p {threads} -G {input.gff} -A {output.genes} -o {output.gtf} -eb {params.outfolder} && "
            "gzip -k {output.gtf} 2> {log}"

if check('fq'):
    rule quantify_transcripts_fq:
        input:
            bam="{dir}/align/{fq}.fq.sorted.bam",
            gff="{dir}/isoforms/combined.gtf",
        output:
            ["{dir}/isoforms/{fq}_quant/t_data.ctab",
            "{dir}/isoforms/{fq}_quant/e_data.ctab",
            "{dir}/isoforms/{fq}_quant/i_data.ctab",
            "{dir}/isoforms/{fq}_quant/e2t.ctab",
            "{dir}/isoforms/{fq}_quant/i2t.ctab"],
            genes="{dir}/isoforms/{fq}_quant/{fq}.fq.gene.quant.tab",
            gtf=temp("{dir}/isoforms/{fq}_quant/{fq}.fq.transcript.quant.gtf"),
            gtfgz="{dir}/isoforms/{fq}_quant/{fq}.fq.transcript.quant.gtf.gz"
        threads: 4
        log: "{dir}/isoforms/{fq}_quant/{fq}.fq.quant.gtf.log"
        params: outfolder=lambda w, output: os.path.dirname(output.genes)
        benchmark: "{dir}/isoforms/{fq}_quant/{fq}.fq.quant.gtf.benchmark"
        conda: "../envs/quant.yaml"
        shell:
            "stringtie {input.bam} -p {threads} -G {input.gff} -A {output.genes} -o {output.gtf} -eb {params.outfolder} && "
            "gzip -k {output.gtf} 2> {log}"

    rule quantify_ref_transcripts_fq:
        input:
            bam="{dir}/align/{fq}.fq.sorted.bam",
            gff=GFF3,
        output:
            ["{dir}/isoforms/{fq}_quant_ref/t_data.ctab",
            "{dir}/isoforms/{fq}_quant_ref/e_data.ctab",
            "{dir}/isoforms/{fq}_quant_ref/i_data.ctab",
            "{dir}/isoforms/{fq}_quant_ref/e2t.ctab",
            "{dir}/isoforms/{fq}_quant_ref/i2t.ctab"],
            genes="{dir}/isoforms/{fq}_quant_ref/{fq}.fq.gene.quant_ref.tab",
            gtf=temp("{dir}/isoforms/{fq}_quant_ref/{fq}.fq.transcript.quant_ref.gtf"),
            gtfgz="{dir}/isoforms/{fq}_quant_ref/{fq}.fq.transcript.quant_ref.gtf.gz"
        threads: 4
        log: "{dir}/isoforms/{fq}_quant_ref/{fq}.fq.quant_ref.gtf.log"
        params: outfolder=lambda w, output: os.path.dirname(output.genes)
        benchmark: "{dir}/isoforms/{fq}_quant_ref/{fq}.fq.quant_ref.gtf.benchmark"
        conda: "../envs/quant.yaml"
        shell:
            "stringtie {input.bam} -p {threads} -G {input.gff} -A {output.genes} -o {output.gtf} -eb {params.outfolder}  && "
            "gzip -k {output.gtf} 2> {log}"

if check('fq_se'):
    rule quantify_transcripts_fq_se:
        input:
            bam="{dir}/align/{fq_se}.fq_se.sorted.bam",
            gff="{dir}/isoforms/combined.gtf",
        output:
            ["{dir}/isoforms/{fq_se}_quant/t_data.ctab",
            "{dir}/isoforms/{fq_se}_quant/e_data.ctab",
            "{dir}/isoforms/{fq_se}_quant/i_data.ctab",
            "{dir}/isoforms/{fq_se}_quant/e2t.ctab",
            "{dir}/isoforms/{fq_se}_quant/i2t.ctab"],
            genes="{dir}/isoforms/{fq_se}_quant/{fq_se}.fq_se.gene.quant.tab",
            gtf=temp("{dir}/isoforms/{fq_se}_quant/{fq_se}.fq_se.transcript.quant.gtf"),
            gtfgz="{dir}/isoforms/{fq_se}_quant/{fq_se}.fq_se.transcript.quant.gtf.gz"
        threads: 4
        log: "{dir}/isoforms/{fq_se}_quant/{fq_se}.fq_se.quant.gtf.log"
        params: outfolder=lambda w, output: os.path.dirname(output.genes)
        benchmark: "{dir}/isoforms/{fq_se}_quant/{fq_se}.fq_se.quant.gtf.benchmark"
        conda: "../envs/quant.yaml"
        shell:
            "stringtie {input.bam} -p {threads} -G {input.gff} -A {output.genes} -o {output.gtf} -eb {params.outfolder} && "
            "gzip -k {output.gtf} 2> {log}"

    rule quantify_ref_transcripts_fq_se:
        input:
            bam="{dir}/align/{fq_se}.fq_se.sorted.bam",
            gff=GFF3,
        output:
            ["{dir}/isoforms/{fq_se}_quant_ref/t_data.ctab",
            "{dir}/isoforms/{fq_se}_quant_ref/e_data.ctab",
            "{dir}/isoforms/{fq_se}_quant_ref/i_data.ctab",
            "{dir}/isoforms/{fq_se}_quant_ref/e2t.ctab",
            "{dir}/isoforms/{fq_se}_quant_ref/i2t.ctab"],
            genes="{dir}/isoforms/{fq_se}_quant_ref/{fq_se}.fq_se.gene.quant_ref.tab",
            gtf=temp("{dir}/isoforms/{fq_se}_quant_ref/{fq_se}.fq_se.transcript.quant_ref.gtf"),
            gtfgz="{dir}/isoforms/{fq_se}_quant_ref/{fq_se}.fq_se.transcript.quant_ref.gtf.gz"
        threads: 4
        log: "{dir}/isoforms/{fq_se}_quant_ref/{fq_se}.fq_se.quant_ref.gtf.log"
        params: outfolder=lambda w, output: os.path.dirname(output.genes)
        benchmark: "{dir}/isoforms/{fq_se}_quant_ref/{fq_se}.fq_se.quant_ref.gtf.benchmark"
        conda: "../envs/quant.yaml"
        shell:
            "stringtie {input.bam} -p {threads} -G {input.gff} -A {output.genes} -o {output.gtf} -eb {params.outfolder} && "
            "gzip -k {output.gtf} 2> {log}"

rule make_isoform_quant_dataframe_custom:
    '''Take the stringtie quant results and put them into a combined dataframe'''
    input:
        lambda w:
            ([] if not check('sra') else expand("{{dir}}/isoforms/{sra}_quant/{sra}.sra.transcript.quant.gtf", sra=config["sra"])) + \
            ([] if not check('sra_se') else expand("{{dir}}/isoforms/{sra_se}_quant/{sra_se}.sra_se.transcript.quant.gtf", sra_se=config["sra_se"])) + \
            ([] if not check('fq') else expand("{{dir}}/isoforms/{fq}_quant/{fq}.fq.transcript.quant.gtf", fq=config["fq"])) + \
            ([] if not check('fq_se') else expand("{{dir}}/isoforms/{fq_se}_quant/{fq_se}.fq_se.transcript.quant.gtf", fq_se=config["fq_se"]))
    output: "{dir}/final/transcript_custom_quant.tpms.csv"
    log: "{dir}/final/transcript_custom_quant.tpms.log"
    benchmark: "{dir}/final/transcript_custom_quant.tpms.benchmark"
    conda: "../envs/quant.yaml"
    shell:
        "python scripts/SummarizeQuantGtf.py {output} {input} &> {log}"

rule make_isoform_quant_dataframe_ref:
    '''Take the stringtie quant results and put them into a combined dataframe'''
    input:
        lambda w:
            ([] if not check('sra') else expand("{{dir}}/isoforms/{sra}_quant_ref/{sra}.sra.transcript.quant_ref.gtf", sra=config["sra"])) + \
            ([] if not check('sra_se') else expand("{{dir}}/isoforms/{sra_se}_quant_ref/{sra_se}.sra_se.transcript.quant_ref.gtf", sra_se=config["sra_se"])) + \
            ([] if not check('fq') else expand("{{dir}}/isoforms/{fq}_quant_ref/{fq}.fq.transcript.quant_ref.gtf", fq=config["fq"])) + \
            ([] if not check('fq_se') else expand("{{dir}}/isoforms/{fq_se}_quant_ref/{fq_se}.fq_se.transcript.quant_ref.gtf", fq_se=config["fq_se"]))
    output: "{dir}/final/transcript_reference_quant.tpms.csv"
    log: "{dir}/final/transcript_reference_quant.tpms.log"
    benchmark: "{dir}/final/transcript_reference_quant.tpms.benchmark"
    conda: "../envs/quant.yaml"
    shell:
        "python scripts/SummarizeQuantGtf.py {output} {input} &> {log}"

rule make_gene_quant_dataframe_custom:
    '''Take the stringtie quant results and put them into a combined dataframe'''
    input:
        lambda w:
            ([] if not check('sra') else expand("{{dir}}/isoforms/{sra}_quant/{sra}.sra.gene.quant.tab", sra=config["sra"])) + \
            ([] if not check('sra_se') else expand("{{dir}}/isoforms/{sra_se}_quant/{sra_se}.sra_se.gene.quant.tab", sra_se=config["sra_se"])) + \
            ([] if not check('fq') else expand("{{dir}}/isoforms/{fq}_quant/{fq}.fq.gene.quant.tab", fq=config["fq"])) + \
            ([] if not check('fq_se') else expand("{{dir}}/isoforms/{fq_se}_quant/{fq_se}.fq_se.gene.quant.tab", fq_se=config["fq_se"]))
    output: "{dir}/final/gene_custom_quant.tpms.csv"
    log: "{dir}/final/gene_custom_quant.tpms.log"
    benchmark: "{dir}/final/gene_custom_quant.tpms.benchmark"
    conda: "../envs/quant.yaml"
    shell:
        "python scripts/SummarizeQuantTab.py {output} {input} &> {log}"

rule make_gene_quant_dataframe_ref:
    '''Take the stringtie quant results and put them into a combined dataframe'''
    input:
        lambda w:
            ([] if not check('sra') else expand("{{dir}}/isoforms/{sra}_quant_ref/{sra}.sra.gene.quant_ref.tab", sra=config["sra"])) + \
            ([] if not check('sra_se') else expand("{{dir}}/isoforms/{sra_se}_quant_ref/{sra_se}.sra_se.gene.quant_ref.tab", sra_se=config["sra_se"])) + \
            ([] if not check('fq') else expand("{{dir}}/isoforms/{fq}_quant_ref/{fq}.fq.gene.quant_ref.tab", fq=config["fq"])) + \
            ([] if not check('fq_se') else expand("{{dir}}/isoforms/{fq_se}_quant_ref/{fq_se}.fq_se.gene.quant_ref.tab", fq_se=config["fq_se"]))
    output: "{dir}/final/gene_reference_quant.tpms.csv"
    log: "{dir}/final/gene_reference_quant.tpms.log"
    benchmark: "{dir}/final/gene_reference_quant.tpms.benchmark"
    conda: "../envs/quant.yaml"
    shell:
        "python scripts/SummarizeQuantTab.py {output} {input} &> {log}"
