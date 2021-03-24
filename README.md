# Spritz
Software for RNA-Seq analysis on Windows, including creating sample-specific proteoform databases from genomic data

[![Build status](https://ci.appveyor.com/api/projects/status/p54yrm6iixqm8jsf?svg=true)](https://ci.appveyor.com/project/acesnik/spritz)
[![Release](https://img.shields.io/github/v/release/smith-chem-wisc/Spritz)](https://github.com/smith-chem-wisc/Spritz/releases/latest)
[![Github All Releases](https://img.shields.io/github/downloads/smith-chem-wisc/Spritz/total.svg)](https://github.com/smith-chem-wisc/Spritz/releases/)
[![Docker Pulls](https://img.shields.io/docker/pulls/smithlab/spritz)](https://hub.docker.com/r/smithlab/spritz/tags?page=1&ordering=last_updated)
[![Follow us on Twitter](https://img.shields.io/twitter/follow/smith_chem_wisc?label=Twitter&style=social)](https://twitter.com/smith_chem_wisc)

Spritz can be downloaded [here](https://github.com/smith-chem-wisc/Spritz/releases).

Spritz uses snakemake and Docker to install and run commandline tools for Next-Generation Sequencing (NGS) analysis. These tools include an [adapted version of SnpEff](https://github.com/smith-chem-wisc/SnpEff) to annotate sequence variations and create an annotated protein database in XML format. The combinatorics of producing full-length proteoforms from these annotations is written in [mzLib's VariantApplication class](https://github.com/smith-chem-wisc/mzLib/blob/master/Proteomics/Protein/VariantApplication.cs).

![image](https://user-images.githubusercontent.com/16342951/93618988-a3b5be00-f99d-11ea-8be4-063395e24ce1.png)

## Running Spritz with GUI

1. Install [Docker Desktop for Windows](https://hub.docker.com/editions/community/docker-ce-desktop-windows).

2. Allocate resources to Docker. There are two ways to do this, described in the Spritz wiki:
    1. The [recommended method](https://github.com/smith-chem-wisc/Spritz/wiki/Allocating-resources-to-Docker-using-WSL2-(recommended-method)) requires Windows 10 version 2004 and is more robust. Here, we allocate computer resources to Docker like any other program.
    2. The [alternate method](https://github.com/smith-chem-wisc/Spritz/wiki/Allocating-resources-to-Docker-using-Vmmem-(alternate-method)) is available on all Windows versions but is less robust. Here, we allocate computer resources to Docker using a virtual machine that's packaged with Docker.

3. Launch Spritz.

    Step 1: Input SRA accessions and/or add FASTQ files.
    * SRAs are added with the button indicating single-end or paired-end.
    * FASTQ files must end with *_1.fastq if single-end, and paired-end sequences must have the same filename other than ending with *_1.fastq and *_2.fastq.

    Step 2: Create and customize your Spritz workflow.

    Step 3: Run Spritz!

    ![Intro-01](https://user-images.githubusercontent.com/16342951/99159330-088d4c00-26a1-11eb-99b9-0b7cd79467d2.png)

## GUI System Requirements

* Environment:
  * Windows 10 recommended
  * [.NET Core 3.1](https://dotnet.microsoft.com/download/dotnet-core/thank-you/runtime-desktop-3.1.3-windows-x64-installer)
* 24 GB RAM recommended
* The installer ([Spritz.msi](https://github.com/smith-chem-wisc/Spritz/releases)) only works on Windows.

## Running Spritz with commandline

Spritz will also [work on the commandline](https://github.com/smith-chem-wisc/Spritz/wiki/Spritz-commandline-usage) within a Unix system (Linux, Mac, WSL on Windows).

## Test it out! Try constructing the database for U2OS from the paper.

1. Add SRR629563 to the SRA list.

2. Create the Spritz workflow. Select "release-97" and "homo_sapiens."

3. Run Spritz!

Monitor progress in the Information textbox. The final database named `final/combined.spritz.snpeff.protein.withmods.xml.gz` can be used to search MS/MS with [MetaMorpheus](https://github.com/smith-chem-wisc/MetaMorpheus) to find variant peptides and proteoforms, possibly with modifications. We recommend performing 1) Calibration, 2) Global PTM Discovery (G-PTM-D), and 3) Search tasks to get the best results.

![image](https://user-images.githubusercontent.com/16342951/85874687-a76be700-b798-11ea-9bff-9f68646b03de.png)

The final database named `final/combined.spritz.snpeff.protein.fasta` is generated to contain variant protein sequences, and it may be used in other search software, such as Proteome Discoverer, ProSight, and MASH Explorer.

The final database named `final/combined.spritz.snpeff.protein.withdecoys.fasta` is ready for use in MSFragger. It is generated to contain variant protein sequences with decoy protein sequences appended.

## Citations

If you use this Spritz, please cite:
* `Spritz`: Cesnik, A. J.; Miller, R. M.; Ibrahim, K.; Lu, L.; Millikin, R. J.; Shortreed, M. R.; Frey, B. L.; Smith, L. M. “Spritz: A Proteogenomic Database Engine.” J. Proteome Res. 2020, in press. https://pubs.acs.org/doi/abs/10.1021/acs.jproteome.0c00407

This pipeline uses the following tools:
  * `sra-toolkit`: Leinonen, R.; et al. International Nucleotide Sequence Database Collaboration. The Sequence Read Archive. Nucleic Acids Res. 2011, 39 (Database issue), D19-21. https://doi.org/10.1093/nar/gkq1019.
  * `fastp`: Chen, S.; et al. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics 2018, 34 (17), i884-i890. https://academic.oup.com/bioinformatics/article/34/17/i884/5093234
  * `hisat2`: Kim, D.; et al. Graph-Based Genome Alignment and Genotyping with HISAT2 and HISAT-Genotype. Nat. Biotechnol. 2019, 37 (8), 907-915. https://doi.org/10.1038/s41587-019-0201-4."
  * `samtools`: Li, H.; et al. The Sequence Alignment/Map format and SAMtools. Bioinformatics 2009, 25 (16), 2078-2079. https://academic.oup.com/bioinformatics/article/25/16/2078/204688.
  * `GATK`: McKenna, A.; et al. The Genome Analysis Toolkit: A MapReduce Framework for Analyzing next-Generation DNA Sequencing Data. Genome Res. 2010, 20 (9), 1297-1303. https://doi.org/10.1101/gr.107524.110.
  * `SnpEff`: Cingolani, P.; et al. A Program for Annotating and Predicting the Effects of Single Nucleotide Polymorphisms, SnpEff: SNPs in the Genome of Drosophila Melanogaster Strain W1118; Iso-2; Iso-3. Fly (Austin) 2012, 6 (2), 80-92. https://doi.org/10.4161/fly.19695.
  * `StringTie2`: Kovaka, S.; et al. Transcriptome assembly from long-read RNA-seq alignments with StringTie2. Genome Biol 2019, 20 (278), 1-13. https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1910-1
