sample_id, = glob_wildcards(config["input_dir"] + "/bc/{sample}.R1.fastq.gz")
count_samples, = glob_wildcards(config["input_dir"] + "/count/{sample}.fastq.gz")

rule all:
  input: config["experiment"] + ".counts.by.variant.txt"

rule merge_fastq:
  input: f=config["input_dir"] + "/bc/{sample}.R1.fastq.gz",
         r=config["input_dir"] + "/bc/{sample}.R2.fastq.gz"
  output: "{sample}.merged.fastq.gz"
  shell:
    """
    module add bbmap
    bbmerge.sh in={input.f} in2={input.r} out={output} minoverlap=75
    """

rule cutadapt_bc:
  input: "{sample}.merged.fastq.gz"
  output: "{sample}.bc.fastq"
  params: adapter_seq=config["adapter_seq"] 
  shell:
    """
    module add cutadapt
    cutadapt -j 24 \
    -a {params.adapter_seq} \
    -O 27 \
    -o {output} {input}
    """

rule cutadapt_oligo:
  input: "{sample}.merged.fastq.gz"
  output: "{sample}.oligo.fastq"
  params: adapter_seq=config["adapter_seq"]
  shell:
    """
    module add cutadapt
    cutadapt -j 24 \
    -g {params.adapter_seq} \
    -O 27 \
    -o {output} {input}
    """

rule generate_design:
  input: config["design_file"]
  output: config["experiment"] + "_design.fa"
  script: "scripts/generate_design.R"

rule map_oligos:
  input: config["experiment"] + "_design.fa",
         "{sample}.bc.fastq",
         "{sample}.oligo.fastq"
  output: "{sample}.sequences.barcodes.txt"
  shell: "./scripts/map_oligos.sh {input[0]} {input[1]} {input[2]} {wildcards.sample}"

rule join_sequences:
  input: expand("{sample}.sequences.barcodes.txt", sample = sample_id)
  output: config["experiment"] + ".join.txt"
  shell: "./scripts/join_barcoded_sequences.sh {input} {output}"

rule barcode_map:
  input: config["experiment"] + ".join.txt"
  output: config["experiment"] + ".barcodes.fasta"
  script: "scripts/generate_barcode_map.R"

rule map_dna_rna:
  input: config["input_dir"] + "/count/{sample}.fastq.gz",
         config["experiment"] + ".barcodes.fasta"
  output: "{sample}.reads.txt"
  shell: "./scripts/match_barcodes.sh {input[0]} {wildcards.sample} {input[1]}"

rule compile_counts:
  input: expand("{sample}.reads.txt", sample = count_samples)
  output: config["experiment"] + ".counts.by.variant.txt"
  script: "scripts/compile_counts.R"
