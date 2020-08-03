import os
import sys
import pysam
from subprocess import CalledProcessError

import pandas as pd

SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)

shell.executable("/bin/bash")
pythonpath_string = "%s/../scripts" % SNAKEMAKE_DIR
shell.prefix("source %s/env.cfg; set -euo pipefail; export PYTHONPATH=%s; " % (SNAKEMAKE_DIR, pythonpath_string))

if config == {}:
	configfile: "mrsfast_config.yaml"

compress_dict = {'cram' : 'rc', 'bam' : 'rb'}


MANIFEST = config["manifest"]
REFERENCE = config["reference"]
MASKED_REF = config[REFERENCE]["masked_ref"]
CONTIGS_FILE = config[REFERENCE]["contigs"]
ALIGNED_REF = config["aligned_reference"]
INPUT_TYPE = config["input_format"]
BAM_PARTITIONS = config["bam_partitions"]
UNMAPPED_PARTITIONS = config["unmapped_partitions"]
if UNMAPPED_PARTITIONS == -1:
	UNMAPPED_PARTITIONS = max(BAM_PARTITIONS // 100, 1)


COMPRESSION = compress_dict[INPUT_TYPE]

MAX_EDIST = config["max_edist"]

TMPDIR = config["tmpdir"]
CLEAN_TEMP_FILES = config["clean_temp_files"]

if not os.path.exists("log"):
	os.makedirs("log")

CONTIGS = {}

with open(CONTIGS_FILE, "r") as reader:
	for line in reader:
		contig, size = line.rstrip().split()[0:2]
		CONTIGS[contig] = int(size)

SAMPLES = pd.read_csv(MANIFEST, sep='\t')
SAMPLES.index = SAMPLES.sn

localrules: all, clean, make_chunker

rule all:
	input:  expand("mapping/{sample}/{sample}/wssd_out_file", sample = SAMPLES.sn)

rule clean:
	input: "mapping/{sample}/{sample}/wssd_out_file"
	output: touch("finished/{sample}.txt")
	priority: 50
	run:
		if CLEAN_TEMP_FILES:
			shell("rm region_matrices/{wildcards.sample}/*")

if config["mode"] == "full":
	rule full_index:
		input: 
			bam = lambda wildcards: SAMPLES.ix[SAMPLES.sn == wildcards.sample, "path"]
		output: 
			index = temp('{sample}_reads.txt')
		params: sge_opts = "-pe serial 1 -l mfree=8G -N map_count -l h_rt=10:00:00"
		resources:
			mem = 16
		threads: 1
		priority: 20
		shell:
			'''
			samtools idxstats {input.bam} | awk '{{sum+=($3+$4)}} END {{print sum}}' > {output.index}
			'''

	rule full_map:
		input: 
			bam = lambda wildcards: SAMPLES.ix[SAMPLES.sn == wildcards.sample, "path"], 
			index = rules.full_index.output.index
		output: 
			tab = "mapped/{sample}/{sample}.{part}_%d.tab.gz" % (BAM_PARTITIONS)
		params: 
			sge_opts = "-pe serial 4 -l mfree=4G -N map_count -l h_rt=10:00:00 -l disk_free=2G"
		benchmark: "benchmarks/counter/{sample}/{sample}.{part}.%d.txt" % BAM_PARTITIONS
		resources: 
			mem=8
		threads: 4
		priority: 20
		log: "log/map/{sample}/{part}_%s.txt" % BAM_PARTITIONS
		run:
			fifo = "$TMPDIR/mrsfast_fifo"
			if TMPDIR != "":
				local_index = "%s/%s" % (TMPDIR, os.path.basename(input.index[0]))
			else:
				local_index = input.index[0]
			shell("hostname; echo part: {wildcards.part} nparts: {BAM_PARTITIONS} unmapped parts: {UNMAPPED_PARTITIONS}; mkfifo {fifo}")
			shell("""python {SNAKEMAKE_DIR}/scripts/align_chunk.py -r {ALIGNED_REF} -f {input.bam} -i {input.index} -s {wildcards.part} -p {BAM_PARTITIONS} -c {COMPRESSION} -o $TMPDIR/{wildcards.sample}.{wildcards.part}.bam;""")
			shell("""
					samtools index $TMPDIR/{wildcards.sample}.{wildcards.part}.bam;
					samtools view -h $TMPDIR/{wildcards.sample}.{wildcards.part}.bam | \
					python {SNAKEMAKE_DIR}/scripts/alignment_to_fasta.py -p {COMPRESSION} -i /dev/stdin -c 36 -o /dev/stdout | \
					mrsfast --threads {threads} --search {MASKED_REF} -n 0 -e {MAX_EDIST} --crop 36 --seq /dev/stdin -o {fifo} --disable-nohit >> /dev/stderr | \
					python {SNAKEMAKE_DIR}/scripts/mrsfast_outputconverter.py {fifo} {output.tab} --compress"""
					)

	rule full_count:
		input: expand("mapped/{{sample}}/{{sample}}.{part}_%d.tab.gz" % (BAM_PARTITIONS), part=range(1, BAM_PARTITIONS + UNMAPPED_PARTITIONS))
		output: "mapping/{sample}/{sample}/wssd_out_file"
		priority: 40
		params: 
			sge_opts="-l mfree=20G -l h_rt=24:00:00 -N map_{sample}"
		shell:
			'''
			python3 {SNAKEMAKE_DIR}/scripts/read_counter_from_gzfile.py {output} --infiles {input} --contigs_file {CONTIGS_FILE}
			'''





elif config["mode"] == "fast":
	rule fast_map:
		input: 
			reads=lambda wildcards: SAMPLES.ix[SAMPLES.sn == wildcards.sample, "path"]
		output: 
			mapped = "mapped/{sample}.tab.gz"
		params: 
			sge_opts="-l mfree=4G -pe serial 16 -l h_rt=4:0:00:00 -N map_{sample}",
			fifo="$TMPDIR/{sample}.fifo"
		benchmark: "benchmarks/{sample}_mapping.txt"
		resources:
			mem = 4
		threads: 16
		shell:
			"""mkfifo {params.fifo};
			   pv -L 50M {input.reads} -q |
			   samtools view -h -T {ALIGNED_REF} - | \
			   python {SNAKEMAKE_DIR}/scripts/sam_to_fastq.py /dev/stdin /dev/stdout --min_length 36 --offset 0 | \
			   mrsfast --search {MASKED_REF} --crop 36 -n 0 -e 2 --seq /dev/stdin -o /dev/stdout \
					   --disable-nohit --threads {threads} --mem 56 | 
			   python {SNAKEMAKE_DIR}/scripts/mrsfast_outputconverter.py /dev/stdin {output.mapped} --compress"""

	rule fast_count:
		input: rules.fast_map.output.mapped
		output: "mapping/{sample}/{sample}/wssd_out_file"
		params: 
			sge_opts="-l mfree=20G -l h_rt=24:00:00 -N map_{sample} -l disk_free=20G"
		resources:
			mem = 20
		threads: 1
		shell:
			"""python3 {SNAKEMAKE_DIR}/scripts/read_counter_from_gzfile.py {output} --infiles {input} --contigs_file {CONTIGS_FILE}"""

else:
	raise ValueError("config['mode'] must be in ['fast', 'full'] but it is {}".format(config.get("mode", "not set")))

rule check_bam_files:
	input: [bam for bam in SAMPLES.path]
	output: touch("BAMS_READABLE")
	params: sge_opts = "-l h_rt=1:0:0"
	priority: 50
	run:
		for bamfile in input:
			try:
				test = pysam.AlignmentFile(bamfile)
			except ValueError as e:
				print("Error: could not open %s as bam.\n%s\n" % (bamfile, str(e)), file = sys.stderr)
				sys.exit(1)

rule check_index:
	input: MASKED_REF
	output: touch("MRSFASTULTRA_INDEXED"), temp(".mrsfast_index_test_output.txt")
	params: sge_opts = "-l h_rt=1:0:0"
	run:
		try:
			shell("mrsfast --search {input[0]} --seq dummy.fq > {output[1]}")
		except CalledProcessError as e:
			sys.exit("Reference %s was not indexed with the current version of mrsfastULTRA. Please reindex." % input[0])

rule make_chunker:
	input: "src/chunker_cascade.cpp", "Makefile"
	output: "bin/bam_chunker_cascade"
	params: sge_opts = "-l h_rt=1:0:0"
	shell:
		"make"

