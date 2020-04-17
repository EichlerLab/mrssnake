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

localrules: all, get_headers, make_jobfile, clean, make_chunker

rule all:
    input:  expand("mapping/{sample}/{sample}/wssd_out_file", sample = SAMPLES.sn)

rule clean:
    input: "mapping/{sample}/{sample}/wssd_out_file"
    output: touch("finished/{sample}.txt")
    priority: 50
    run:
        if CLEAN_TEMP_FILES:
            shell("rm region_matrices/{wildcards.sample}/*")

if config["mode"] == "full" and INPUT_TYPE == 'bam':
    rule full_count:
        input: expand("mapped/{{sample}}/{{sample}}.{part}_%d.tab.gz" % (BAM_PARTITIONS), part=range(BAM_PARTITIONS + UNMAPPED_PARTITIONS))
        output: "mapping/{sample}/{sample}/wssd_out_file"
        params: sge_opts="-l mfree=20G -l h_rt=24:00:00 -N map_{sample} -l disk_free=20G"
        shell:
            '''
            python3 scripts/read_counter_from_gzfile.py {output} --infiles {input} --contigs_file {CONTIGS_FILE}
            '''

    rule full_map:
        input: bam = lambda wildcards: SAMPLES.ix[SAMPLES.sn == wildcards.sample, "path"], index = lambda wildcards: SAMPLES.ix[SAMPLES.sn == wildcards.sample, "index"]
        output: "mapped/{sample}/{sample}.{part}_%d.tab.gz" % (BAM_PARTITIONS)
        params: sge_opts = "-pe serial 4 -l mfree=4G -N map_count -l h_rt=10:00:00"
        benchmark: "benchmarks/counter/{sample}/{sample}.{part}.%d.txt" % BAM_PARTITIONS
        resources: mem=10
        priority: 20
        log: "log/map/{sample}/{part}_%s.txt" % BAM_PARTITIONS
        run:
            masked_ref_name = os.path.basename(MASKED_REF)
            fifo = "%s/mrsfast_fifo" % TMPDIR
            if TMPDIR != "":
                local_index = "%s/%s" % (TMPDIR, os.path.basename(input.index[0]))
            else:
                local_index = input.index[0]
            mrsfast_ref_path = "/var/tmp/$USER/%s" % masked_ref_name
            rsync_opts = """mkdir -p /var/tmp/$USER; 
                            rsync {0}.index /var/tmp/$USER/ --bwlimit 10000 --copy-links -p;
                            rsync {2} {3} --bwlimit 10000 --copy-links -p; 
                            if [[ ! -e {1} ]]; then touch {1}; fi; 
                            echo Finished rsync from {0} to {1} >> /dev/stderr; 
                            echo Finished rsync from {2} to {3} >> /dev/stderr; """.format(MASKED_REF, mrsfast_ref_path, input.index[0], local_index)

            shell("hostname; echo part: {wildcards.part} nparts: {BAM_PARTITIONS} unmapped parts: {UNMAPPED_PARTITIONS}; mkfifo {fifo}; {rsync_opts}")
            shell("{SNAKEMAKE_DIR}/bin/bam_chunker -b {input.bam} -i {local_index} -p {wildcards.part} -n {BAM_PARTITIONS} -u {UNMAPPED_PARTITIONS} 2>> /dev/stderr | "
                "mrsfast --search {mrsfast_ref_path} -n 0 -e {MAX_EDIST} --crop 36 --seq /dev/stdin -o {fifo} --disable-nohit >> /dev/stderr | "
                "python3 scripts/mrsfast_outputconverter.py {fifo} {output} --compress"
                )

elif config["mode"] == "full" and INPUT_TYPE == 'cram':
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

    rule full_map:
        input: bam = lambda wildcards: SAMPLES.ix[SAMPLES.sn == wildcards.sample, "path"], index = '{sample}_dump.txt'
        output: tab = "mapped/{sample}/{sample}.{part}_%d.tab.gz" % (BAM_PARTITIONS)
        params: sge_opts = "-pe serial 4 -l mfree=4G -N map_count -l h_rt=10:00:00 -l disk_free=2G"
        benchmark: "benchmarks/counter/{sample}/{sample}.{part}.%d.txt" % BAM_PARTITIONS
        resources: mem=10
        priority: 20
        log: "log/map/{sample}/{part}_%s.txt" % BAM_PARTITIONS
        run:
            fifo = "$TMPDIR/mrsfast_fifo"
            if TMPDIR != "":
                local_index = "%s/%s" % (TMPDIR, os.path.basename(input.index[0]))
            else:
                local_index = input.index[0]
            shell("hostname; echo part: {wildcards.part} nparts: {BAM_PARTITIONS} unmapped parts: {UNMAPPED_PARTITIONS}; mkfifo {fifo}")
            shell("""python {SNAKEMAKE_DIR}/scripts/cramChunker.py -c {input.bam} -i {input.index} -s {wildcards.part} -p {BAM_PARTITIONS} -o $TMPDIR/{wildcards.sample}.{wildcards.part}.cram; \
                samtools index $TMPDIR/{wildcards.sample}.{wildcards.part}.cram;
                samtools view -h -T {ALIGNED_REF} $TMPDIR/{wildcards.sample}.{wildcards.part}.cram | \
                python {SNAKEMAKE_DIR}/scripts/cram_to_fasta.py -i /dev/stdin -c 36 -o /dev/stdout | \
                mrsfast --search {MASKED_REF} -n 0 -e {MAX_EDIST} --crop 36 --seq /dev/stdin -o {fifo} --disable-nohit >> /dev/stderr | \
                python3 {SNAKEMAKE_DIR}/scripts/mrsfast_outputconverter.py {fifo} {output.tab} --compress"""
                )

    rule full_index:
        input: 
            bam = lambda wildcards: SAMPLES.ix[SAMPLES.sn == wildcards.sample, "path"]
        output: 
            index = temp('{sample}_dump.txt')
        params: sge_opts = "-pe serial 1 -l mfree=8G -N map_count -l h_rt=10:00:00"
        priority: 20
        shell:
            '''
            cram_dump {input.bam} | grep "Container pos" > {output.index}
            '''


elif config["mode"] == "fast":
    rule fast_count:
        input: "mapped/{sample}.tab.gz"
        output: "mapping/{sample}/{sample}/wssd_out_file"
        params: sge_opts="-l mfree=20G -l h_rt=24:00:00 -N map_{sample} -l disk_free=20G"
        shell:
            """python3 {SNAKEMAKE_DIR}/scripts/read_counter_from_gzfile.py {output} --infiles {input} --contigs_file {CONTIGS_FILE}"""

    rule fast_map:
        input: reads=lambda wildcards: SAMPLES.ix[SAMPLES.sn == wildcards.sample, "path"]
        output: "mapped/{sample}.tab.gz"
        params: sge_opts="-l mfree=4G -pe serial 4 -l h_rt=4:0:00:00 -N map_{sample}",
                fifo="$TMPDIR/{sample}.fifo"
        benchmark: "benchmarks/{sample}_mapping.txt"
        shell:
            """mkfifo {params.fifo};
               pv -L 50M {input.reads} -q |
               samtools view -h -T {ALIGNED_REF} - | \
               python {SNAKEMAKE_DIR}/scripts/sam_to_fastq.py /dev/stdin /dev/stdout --min_length 36 --offset 0 | \
               mrsfast --search {MASKED_REF} --crop 36 -n 0 -e 2 --seq /dev/stdin -o /dev/stdout \
                       --disable-nohit --threads 4 --mem 8 | 
               python {SNAKEMAKE_DIR}/scripts/mrsfast_outputconverter.py /dev/stdin {output} --compress"""
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
