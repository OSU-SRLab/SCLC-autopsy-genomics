import os
import pathlib

cwd=os.getcwd()

indir = pathlib.Path("input")
paths = indir.glob("*.fastq")
#study, patient, T/N, DNA ID
samples = set(["-".join(x.stem.split("-")[0:4]) for x in paths])
gz_paths = indir.glob("*.fastq.gz")
[samples.add("-".join(x.stem.split("-")[0:4])) for x in gz_paths]

trimmed_fastqs_path = "output/fastq_trim"

### substitute your path here
workdir: "/parh/to/analysis_dir"

def find_fastqs(wildcards, read, custom_path = None):
    raw_fastq_gzs = [str(x) for x in indir.glob(f"{wildcards.sample}*_R{read}_*.fastq*")]
    raw_fastqs = set()
    #strip .gz extension to couple properly with rule gunzip_fastq
    for raw_fastq_gz in raw_fastq_gzs:
        if raw_fastq_gz[-3:] == ".gz":
            raw_fastqs.add(raw_fastq_gz[:-3])
        else:
            raw_fastqs.add(raw_fastq_gz)
    fastqs = None
    if custom_path is not None:
        fastqs = []
        for path in raw_fastqs:
            base = os.path.basename(path)
            new_path = os.path.join(custom_path, base)
            fastqs.append(new_path)
    else:
        fastqs = [fastq for fastq in raw_fastqs]
    return sorted(fastqs)

def find_fastqs_1(wildcards):
    return find_fastqs(wildcards, "1")

def find_fastqs_2(wildcards):
    return find_fastqs(wildcards, "2")

#these return paths in the trimmed fastq directory
def find_trimmed_fastqs_1(wildcards):
    return find_fastqs(wildcards, "1", custom_path=trimmed_fastqs_path)

def find_trimmed_fastqs_2(wildcards):
    return find_fastqs(wildcards, "2", custom_path=trimmed_fastqs_path)

def mem_per_thread(wildcards, input, threads):
    return 1024*4*threads

def procs_per_thread(wildcards, input, threads):
    return threads

rule all:
    input:
        expand("output/stringtie_merged_transcript/{sample}_out.gtf", sample=samples),
#        expand("output/bam/{sample}.bam", sample=samples),

#output of this marked as temp since no need to retain both gzip and uncompressed
rule gunzip_fastq:
    input:
        "input/{gzfile}.fastq.gz"
    output:
        temp("input/{gzfile}.fastq")
    threads: 1
    resources:
        runtime="00:45:00",
        mem_mb=mem_per_thread,
        cpus=procs_per_thread,
        nodes=1
    shell:
        "gunzip -c {input} > {output}"

#also does fastqc for original fastqs since no good way to set up dependency to rule all
rule trim_galore:
    input:
        fastq_1 = find_fastqs_1,
        fastq_2 = find_fastqs_2
    output:
        touch("output/flags/trim_galore/{sample}.done")
    params:
        fastqc_orig_output_dir = "output/fastqc",
        rename_fastq_ext = "resources/software/rename_fastq_ext.py"
    threads: 4
    resources:
        runtime="4:00:00",
        mem_mb=mem_per_thread,
        cpus=procs_per_thread,
        nodes=1
    shell:
        "fastqc -o {params.fastqc_orig_output_dir} {input.fastq_1} ; "
        "fastqc -o {params.fastqc_orig_output_dir} {input.fastq_2} ; "
        "trim_galore -j {threads} --fastqc --paired {input.fastq_1} {input.fastq_2} --output_dir " + trimmed_fastqs_path + " ; "
        "python {params.rename_fastq_ext} {input.fastq_1} {input.fastq_2} " + trimmed_fastqs_path

rule hisat2:
    input:
        "output/flags/trim_galore/{sample}.done"
    output:
        mapped = "output/bam/{sample}.bam",
        mapped_idx = "output/bam/{sample}.bam.bai"
    params:
        genome_basename = "genome_snp_tran",
        fastq_1 = find_trimmed_fastqs_1,
        fastq_2 = find_trimmed_fastqs_2,
        sam = "temp/{sample}.sam"
    threads: 6
    resources:
        runtime="7:00:00",
        mem_mb=1024*4*8,
        cpus=8,
        nodes=1
    shell:
        "cd " + cwd + " ; "
        "export HISAT2_INDEXES=" + os.path.join(cwd, "resources/reference_data/hisat2/grch37_snp_tran") + " ; "
        "hisat2 --dta -t -p {threads} -x {params.genome_basename} -1 {params.fastq_1} -2 {params.fastq_2} -S {params.sam} ; "
        "samtools sort -o {output.mapped} {params.sam} ; "
        "rm {params.sam} ; "
        "samtools index {output.mapped}"

rule filter_bam:
    input:
        mapped = "output/bam/{sample}.bam"
    output:
        mapped_filtered = "output/bam_filtered/{sample}.bam",
        mapped_filtered_idx = "output/bam_filtered/{sample}.bam.bai"
    params:
        filter_bed_file = "resources/reference_data/filter_bam/IDT_xGen_Exome_Research_CNV-Core_merged_sorted_grch37.bed",
        genome_lengths_file = "resources/reference_data/filter_bam/grch37_lengths.genome"
    threads: 1
    resources:
        runtime="12:00:00",
        mem_mb=1024*4*3,
        cpus=3,
        nodes=1
    shell:
        "bedtools intersect -a {input.mapped} -b {params.filter_bed_file} -wa -sorted -g {params.genome_lengths_file} > {output.mapped_filtered} ; "
        "samtools index {output.mapped_filtered}"

rule stringtie:
    input:
        mapped_filtered = "output/bam_filtered/{sample}.bam"
    output:
        stringtie_out = "output/stringtie/{sample}_out.gtf",
        stringtie_abund = "output/stringtie/{sample}_abund.tab",
        stringtie_cov_refs = "output/stringtie/{sample}_cov_refs.gtf"
    params:
        stringtie_ref_transcripts = "resources/reference_data/stringtie/Homo_sapiens.GRCh37.87.gtf"
    threads: 3
    resources:
        runtime="2:00:00",
        mem_mb=mem_per_thread,
        cpus=procs_per_thread,
        nodes=1
    shell:
        "stringtie {input.mapped_filtered} -p {threads} -o {output.stringtie_out} -G {params.stringtie_ref_transcripts} -A {output.stringtie_abund} -C {output.stringtie_cov_refs}"

rule stringtie_merge_transcripts:
    input:
        expand("output/stringtie/{sample}_out.gtf", sample=samples)
    output:
        merged_transcripts_gtf = "output/stringtie_merged_transcript/merged_transcripts.gtf"
    params:
        stringtie_ref_transcripts = "resources/reference_data/stringtie/Homo_sapiens.GRCh37.87.gtf"
    threads: 6
    resources:
        runtime="8:00:00",
        mem_mb=mem_per_thread,
        cpus=procs_per_thread,
        nodes=1
    shell:
        "stringtie --merge -p {threads} -G {params.stringtie_ref_transcripts} -o {output.merged_transcripts_gtf} {input}"

rule gffcompare:
    input:
        merged_transcripts_gtf = "output/stringtie_merged_transcript/merged_transcripts.gtf"
    output:
        gffcmp_anno = "output/stringtie_merged_transcript/gffcmp.annotated.gtf"
    params:
        stringtie_ref_transcripts = "resources/reference_data/stringtie/Homo_sapiens.GRCh37.87.gtf"
    threads: 1
    resources:
        runtime="4:00:00",
        mem_mb=mem_per_thread,
        cpus=procs_per_thread,
        nodes=1
    shell:
        "gffcompare -r {params.stringtie_ref_transcripts} -o output/stringtie_merged_transcript/gffcmp {input.merged_transcripts_gtf}"

rule stringtie_with_merged_transcripts:
    input:
        mapped_filtered = "output/bam_filtered/{sample}.bam",
        merged_transcripts_gtf = "output/stringtie_merged_transcript/gffcmp.annotated.gtf"
    output:
        stringtie_out_merged = "output/stringtie_merged_transcript/{sample}_out.gtf",
        stringtie_abund_merged = "output/stringtie_merged_transcript/{sample}_abund.tab",
        stringtie_cov_refs_merged = "output/stringtie_merged_transcript/{sample}_cov_refs.gtf"
    threads: 3
    resources:
        runtime="2:00:00",
        mem_mb=mem_per_thread,
        cpus=procs_per_thread,
        nodes=1
    shell:
        "stringtie {input.mapped_filtered} -p {threads} -o {output.stringtie_out_merged} -G {input.merged_transcripts_gtf} -A {output.stringtie_abund_merged} -C {output.stringtie_cov_refs_merged} -e"
