[![Snakemake](https://img.shields.io/badge/snakemake-≥3.5.4-brightgreen.svg?style=flat-square)](http://snakemake.bitbucket.org)

# mrssnake
Simple mrsfast read depth mapping using snakemake

## Quick start
1. Clone the repository, set up the config file, compile the C extension

   ```bash
   git clone --recursive https://github.com/eichlerlab/mrssnake
   cd mrssnake/scripts
   ./setup.sh
   cd ..
   ```
2. Create a tab-delimited manifest file with the appropriate header and a line for each sample

    sn  | path | index 
    --- | ---- | ----- 
    sample_name  | /full/path/to/reads | /full/path/to/reads_index 


3. Modify `mrssfast_config.yaml`

   In particular, set the manifest variable to point to your manifest file, 
   make sure the `reference` and `aligned_reference` variables point to the appropriate reference, 
   and that the paths for the `masked_ref` and `contigs` for your reference are correct.
   
   To add a new reference `new_ref` to the config file, change `reference`:
   ```yaml
   reference: new_ref
   ```
   And add lines for `new_ref`:
   ```yaml
   new_ref:
       masked_ref: /path/to/masked_ref
       contigs: /path/to/contigs
   ```
   Note that every `masked_ref` must be TRF and repeat-masked and indexed with mrsfastULTRA.
   
4. Run snakemake

   This example will use 100 cores with drmaa:
   ```bash
   export DRMAA_LIBRARY_PATH=/opt/uge/lib/lx-amd64/libdrmaa.so.1.0
   snakemake --drmaa " -V -cwd -e ./log -o ./log {params.sge_opts} -w n -S /bin/bash" -w 30 -j 100 -kT
   ```


