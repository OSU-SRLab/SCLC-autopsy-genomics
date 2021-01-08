# SCLC-autopsy-genomics
Scripts associated with our study, "Genomic and transcriptomic characterization of relapsed small cell lung cancer through rapid research autopsy"

### Dependencies

|--------------|----------------------------|
| **Software** | **URL** | **Version Used** |
|--------------|---------|------------------|
| FASTQC       | http://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc | 0.11.9 |
| HISAT2       | https://daehwankimlab.github.io/hisat2/download/ | 2.2.1 |
| TrimGalore   | https://github.com/FelixKrueger/TrimGalore/releases | 0.6.6 |
| Cutadapt     | https://cutadapt.readthedocs.io/en/stable/ | 3.1 |
| StringTie    | https://ccb.jhu.edu/software/stringtie/#install | 2.1.4 |
| GffCompare | https://github.com/gpertea/gffcompare/releases/tag/v0.12.1 | 0.12.1 |
| SAMtools     | http://www.htslib.org/ | 1.9 |
| bedtools     | https://bedtools.readthedocs.io/en/latest/ | 2.29.2 |
| Snakemake | https://snakemake.readthedocs.io/ | 5.31.1 |
|--------------|---------|------------------|

All executables for the above (`fastqc`, `hisat2`, `trim_galore`, `cutadapt`, `stringtie`, `gffcompare`, `samtools`, `bedtools`) must be available in `$PATH`, and the Python libraries associated with Cutadapt and Snakemake in `$PYTHONPATH`.

### Input file names

This workflow requires paired-end fastqs per sample, named as follows:

```StudyID-PatientNum-NT-SampleID-HybID_Instrument_LaneID_R12.fastq```

where:
|-------|-------------|
| **Field** | **Description** |
|-------|-------------|
| StudyID | String identifying a sequencing study |
| PatientNum | Number or string unique to a patient |
| NT | Either "N" for normal or "T" for tumor |
| SampleID | Number that identifies a tumor sample |
| HybID | Number that identifies a particular hybridization or sequencing run |
| Instrument | Name of sequencing instrument |
| LaneID | Lane number that sample was run on |
| R12 | Either "R1" for read 1 or "R2" for read 2 from paired-end sequencing |
|-------|-------------|

### Running

After cloning the repository, edit the `workdir` in `Snakefile` to point to your working directory. `download_resources.sh` will automatically download reference files for HISAT2 and StringTie. Note that this workflow aligns to hg19 and uses corresponding hg19 reference files; it can be easily adapted to other genome assemblies.

The included BED file (`resources/reference_data/filter_bam/IDT_xGen_Exome_Research_CNV-Core_merged_sorted_grch37.bed`) was derived from the IDT xGen Exome v1 capture region (https://www.idtdna.com/pages/products/next-generation-sequencing/targeted-sequencing/hybridization-capture/predesigned-panels/xgen-exome-research-panel-v2) by removing the "chr" prefix from chromosomes, and merging overlapping regions and sorting with bedtools. Aligned reads from all samples used in our study were filtered with this BED file as we utilized an exome capture in our RNAseq prep, and other studies did not.

This workflow was run at the Ohio Supercomputer Center (https://www.osc.edu/), which utilizes the Slurm batch scheduler (https://slurm.schedmd.com/). We installed the Snakemake slurm profile (https://github.com/Snakemake-Profiles/slurm) using cookiecutter (https://github.com/cookiecutter/cookiecutter). The workflow was then started with Snakemake:
```snakemake --profile slurm --jobs 20```
