rule download_protein_xml:
    '''Download the uniprot xml database and uniprot isoform fasta'''
    output:
        xml=UNIPROTXML,
        fasta=UNIPROTFASTA,
    log: f"{UNIPROTXML}.log"
    benchmark: f"{UNIPROTXML}.benchmark"
    conda: "../envs/downloads.yaml"
    shell:
        "(python scripts/get_proteome.py && "
        "python scripts/download_uniprot.py xml | gzip -c > {output.xml} && " #fixme
        "python scripts/download_uniprot.py fasta > {output.fasta}) &> {log}"

rule build_transfer_mods:
    '''Build the transfer mods C# project'''
    output: TRANSFER_MOD_DLL
    log: "../resources/TransferUniProtModifications.build.log"
    benchmark: "../resources/TransferUniProtModifications.build.benchmark"
    conda: "../envs/proteogenomics.yaml"
    shell:
        "(cd ../TransferUniProtModifications && "
        "dotnet restore && "
        "dotnet build -c Release TransferUniProtModifications.sln) &> {log}"

rule setup_transfer_mods:
    '''Download the ptmlists to the resources directory'''
    input: TRANSFER_MOD_DLL
    output:
        "../resources/ptmlist.txt",
        "../resources/PSI-MOD.obo.xml"
    log: "../resources/setup_transfer_mods.log"
    benchmark: "../resources/setup_transfer_mods.benchmark"
    conda: "../envs/proteogenomics.yaml"
    shell: "cd ../resources/ && dotnet {input} --setup &> {log}"

rule setup_ptmlist_links:
    '''Link the resources to the workflow directory temporarily'''
    input:
        "../resources/ptmlist.txt",
        "../resources/PSI-MOD.obo.xml"
    output:
        temp("ptmlist.txt"),
        temp("PSI-MOD.obo.xml")
    log: "../resources/setup_transfer_mod_linking.log"
    benchmark: "../resources/setup_transfer_mod_linking.benchmark"
    conda: "../envs/proteogenomics.yaml"
    shell: "ln -s {input} . &> {log}"

rule transfer_modifications_variant:
    input:
        "ptmlist.txt", "PSI-MOD.obo.xml",
        transfermods=TRANSFER_MOD_DLL,
        unixml=UNIPROTXML,
        protxml="{dir}/variants/combined.spritz.snpeff.protein.xml",
    output:
        protfastawithdecoys="{dir}/variants/combined.spritz.snpeff.protein.withdecoys.fasta",
        protxmlgz="{dir}/variants/combined.spritz.snpeff.protein.xml.gz",
        protxmlwithmods=temp("{dir}/variants/combined.spritz.snpeff.protein.withmods.xml"),
        protxmlwithmodsgz="{dir}/variants/combined.spritz.snpeff.protein.withmods.xml.gz",
    log: "{dir}/variants/combined.spritz.snpeff.protein.withmods.log"
    benchmark: "{dir}/variants/combined.spritz.snpeff.protein.withmods.benchmark"
    conda: "../envs/proteogenomics.yaml"
    shell:
        "(dotnet {input.transfermods} -x {input.unixml} -y {input.protxml} && "
        "gzip -k {input.protxml} {output.protxmlwithmods}) &> {log}"

rule transfer_modifications_isoformvariant:
    input:
        "ptmlist.txt", "PSI-MOD.obo.xml",
        transfermods=TRANSFER_MOD_DLL,
        unixml=UNIPROTXML,
        protxml="{dir}/variants/combined.spritz.isoformvariants.protein.xml",
    output:
        protfastawithdecoys="{dir}/variants/combined.spritz.isoformvariants.protein.withdecoys.fasta",
        protxmlgz="{dir}/variants/combined.spritz.isoformvariants.protein.xml.gz",
        protxmlwithmods=temp("{dir}/variants/combined.spritz.isoformvariants.protein.withmods.xml"),
        protxmlwithmodsgz="{dir}/variants/combined.spritz.isoformvariants.protein.withmods.xml.gz",
    log: "{dir}/variants/combined.spritz.isoformvariants.protein.withmods.log"
    conda: "../envs/proteogenomics.yaml"
    shell:
        "(dotnet {input.transfermods} -x {input.unixml} -y {input.protxml} && "
        "gzip -k {output.protxmlwithmods} {input.protxml}) &> {log}"

rule generate_reference_snpeff_database:
    input:
        jar="../resources/SnpEff/snpEff.jar",
        gff3=GFF3,
        pfa=f"../resources/ensembl/{REF}.pep.all.fa",
        gfa=KARYOTYPIC_GENOME_FA,
    output:
        pfa=f"../resources/SnpEff/data/{REF}/protein.fa",
        gff3=f"../resources/SnpEff/data/{REF}/genes.gff",
        gfa=f"../resources/SnpEff/data/genomes/{REF}.fa",
        done=f"../resources/SnpEff/data/{REF}/done{REF}.txt",
    resources: mem_mb=16000
    params:
        snpeff_folder=lambda w, input: os.path.dirname(input.jar),
        ref=REF,
        genome_version=GENOME_VERSION
    benchmark: f"../resources/SnpEff/data/{REF}/snpeffdatabase.benchmark"
    log: f"../resources/SnpEff/data/{REF}/snpeffdatabase.log"
    conda: "../envs/proteogenomics.yaml"
    shell:
        "cp {input.gff3} {output.gff3} && "
        "cp {input.pfa} {output.pfa} && "
        "cp {input.gfa} {output.gfa} && "
        "echo \"\n# {params.ref}\" >> {params.snpeff_folder}/snpEff.config && "
        "echo \"{params.ref}.genome : Human genome {params.genome_version} using RefSeq transcripts\" >> {params.snpeff_folder}/snpEff.config && "
        "echo \"{params.ref}.reference : ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/\" >> {params.snpeff_folder}/snpEff.config && "
        "echo \"\t{params.ref}.M.codonTable : Vertebrate_Mitochondrial\" >> {params.snpeff_folder}/snpEff.config && "
        "echo \"\t{params.ref}.MT.codonTable : Vertebrate_Mitochondrial\" >> {params.snpeff_folder}/snpEff.config && "
        "(java -Xmx{resources.mem_mb}M -jar {input.jar} build -gff3 -v {params.ref}) &> {log} && touch {output.done}"

rule reference_protein_xml:
    """
    Create protein XML with sequences from the reference gene model.
    """
    input:
        "ptmlist.txt", "PSI-MOD.obo.xml",
        f"../resources/SnpEff/data/{REF}/done{REF}.txt",
        snpeff="../resources/SnpEff/snpEff.jar",
        fa=f"../resources/ensembl/{REF}.dna.primary_assembly.karyotypic.fa",
        transfermods=TRANSFER_MOD_DLL,
        unixml=UNIPROTXML,
    output:
        done=os.path.join("{dir}/variants/", f"done{REF}.{ENSEMBL_VERSION}.txt"),
        protxml=temp(os.path.join("{dir}/variants/", f"{REF}.{ENSEMBL_VERSION}.protein.xml")),
        protxmlgz=os.path.join("{dir}/variants/", f"{REF}.{ENSEMBL_VERSION}.protein.xml.gz"),
        protfa=os.path.join("{dir}/variants/", f"{REF}.{ENSEMBL_VERSION}.protein.fasta"),
        protwithdecoysfa=os.path.join("{dir}/variants/", f"{REF}.{ENSEMBL_VERSION}.protein.withdecoys.fasta"),
        protxmlwithmods=temp(os.path.join("{dir}/variants/", f"{REF}.{ENSEMBL_VERSION}.protein.withmods.xml")),
        protxmlwithmodsgz=os.path.join("{dir}/variants/", f"{REF}.{ENSEMBL_VERSION}.protein.withmods.xml.gz"),
    params: ref=REF,
    resources: mem_mb=16000
    benchmark: os.path.join("{dir}/variants/", f"{REF}.{ENSEMBL_VERSION}.spritz.benchmark")
    log: os.path.join("{dir}/variants/", f"{REF}.{ENSEMBL_VERSION}.spritz.log")
    conda: "../envs/proteogenomics.yaml"
    shell:
        "(java -Xmx{resources.mem_mb}M -jar {input.snpeff} -v -nostats"
        " -xmlProt {output.protxml} {params.ref} && " # no isoforms, no variants
        "dotnet {input.transfermods} -x {input.unixml} -y {output.protxml} && "
        "gzip -k {output.protxmlwithmods} {output.protxml}) &> {log} && touch {output.done}"

rule custom_protein_xml:
    """
    Create protein XML with sequences from the isoform discovery gene model.
    """
    input:
        "ptmlist.txt", "PSI-MOD.obo.xml",
        snpeff="../resources/SnpEff/snpEff.jar",
        fa=KARYOTYPIC_GENOME_FA,
        isoform_reconstruction=[
            "../resources/SnpEff/data/combined.transcripts.genome.gff3/genes.gff",
            "../resources/SnpEff/data/combined.transcripts.genome.gff3/protein.fa",
            "../resources/SnpEff/data/genomes/combined.transcripts.genome.gff3.fa",
            "../resources/SnpEff/data/combined.transcripts.genome.gff3/done.txt"],
        transfermods=TRANSFER_MOD_DLL,
        unixml=UNIPROTXML,
    output:
        protxml=temp("{dir}/isoforms/combined.spritz.isoform.protein.xml"),
        protwithdecoysfa="{dir}/isoforms/combined.spritz.isoform.protein.withdecoys.fasta",
        protxmlgz="{dir}/isoforms/combined.spritz.isoform.protein.xml.gz",
        protxmlwithmods=temp("{dir}/isoforms/combined.spritz.isoform.protein.withmods.xml"),
        protxmlwithmodsgz="{dir}/isoforms/combined.spritz.isoform.protein.withmods.xml.gz",
        protfa="{dir}/isoforms/combined.spritz.isoform.protein.fasta",
    params:
        ref="combined.transcripts.genome.gff3", # with isoforms
    resources: mem_mb=16000
    benchmark: "{dir}/isoforms/combined.spritz.isoform.benchmark"
    log: "{dir}/isoforms/combined.spritz.isoform.log"
    conda: "../envs/proteogenomics.yaml"
    shell:
        "(java -Xmx{resources.mem_mb}M -jar {input.snpeff} -v -nostats"
        " -xmlProt {output.protxml} {params.ref} < /dev/null && " # isoforms, no variants
        "dotnet {input.transfermods} -x {input.unixml} -y {output.protxml} && "
        "gzip -k {output.protxmlwithmods} {output.protxml}) &> {log}"
