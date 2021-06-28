# Multiple-Sample IntraTumor Heterogeneity Estimator (MITHE)
MITHE implements a somatic variant post-processing pipeline using a system of Perl and Bash scripts. MITHE can be used to call SNV variants and estimate the intratumor heterogeneity between multiple (>2) samples per individual, given a set of user-specified parameters. MITHE can also be used to optimize the parameters that regulate its execution empirically, using technical replicates (same DNA sample sequenced twice independently). MITHE is intended to be run on an HPC environment, interacting directly with a workload manager supporting dependencies between jobs. We developed MITHE on a SLURM environment, but its flexible configuration should make it easy to run it under other workload managers of similar characteristics.

# Installation
We developed MITHE to run on a Linux environment and have only tested it on Linux. Once dependencies are resolved, MITHE's installation only requires downloading the source code and configuring it to run on the user's HPC environment, either manually or using a helper script included in the repository, configAssistant.sh.

## Dependencies
The following programs and libraries need to be available in the HPC system. MITHE provides an easy way of loading modules and executing code before the different dependencies are executed, increasing its flexibility. After the module is loaded and the associated code is run (if necessary), the program executable should be in the environment's PATH. MITHE provides a helper script to configure the system, configAssistant.pl, but more advanced users can configure it manually by editing a configuration file, config.txt.

### Programs
- Annovar
- BEDOPS 
- GATK with UnifiedGenotyper
- GNU parallel
- Platypus
- Samtools (Only if the bam files are not indexed)
- SNPSIFT
- vcftools

#### Other programs that should be available within the Linux environment
- AWK
- BASH
- Perl
- readlink

### Perl libraries:
- Bio::DB::HTS::Tabix 
- Cwd
- Env
- File::Basename
- File::Copy
- Getopt::Long
- Parallel::Loops (Only for fine-grain parallelization)
- Sort::Key::Maker
- Sort::Key::Natural

*NOTE*: If using an HPC system without administrative privileges, these can be installed in a local Perl library location and load them using Mlocal::lib in the Perl configurable executable file. 

## Configuration
### configAssistant.pl
The easiest way of configuring MITHE for a given environment is to run the assistant configuration script. This script will ask the user a series of questions about their system and run some checks on the answers. It also generates a configuration file (config.txt) and places it in MITHE's root directory for its usage.
### Manual preparation of configuration file
Advanced users may want to prepare the configuration file by themselves or edit the one generated by configAssistant.pl to use as a template. For reference, this configuration file is executed in bash using `source`, so that all environment variables are then inherited by the jobs submitted to the workload manager. The environment variables used by MITHE are explained below.
#### Exhaustive list of MITHE Environment variables
##### Required variables
These variables are necessary to indicate the location of different components of the pipeline.
- *MITHE_HOME*: Directory that contains MITHE repository, where MITHE_loop.sh is located.
- *MITHE_INT*: Directory where internal MITHE scripts are located, typically, $MITHE_HOME/internal
- *MITHE_HUMANDB_DIR*: Directory where annovar stores its reference genome
- *MITHE_HUMAN_GENOME*: Human reference genome fasta file
- *MITHE_GNOMAD*: GNOMAD Population allele frequency information [(see gnomAD for instructions on how to generate this file)](#gnomad).

##### Modules
These variables indicate the name of the module that must be loaded before the execution of the program with the same name.
- MITHE_MODULE_PERL
- MITHE_MODULE_ANNOVAR
- MITHE_MODULE_PLATYPUS
- MITHE_MODULE_VCFTOOLS
- MITHE_MODULE_GNUPARALLEL
- MITHE_MODULE_GATK
- MITHE_MODULE_BEDOPS
- MITHE_MODULE_SNPSIFT
- MITHE_MODULE_SAMTOOLS

##### Exes
The string contained on these variables is executed using `eval` before the execution of the program with the same name. To be used if needed.
- MITHE_EXE_PERL
- MITHE_EXE_ANNOVAR
- MITHE_EXE_PLATYPUS
- MITHE_EXE_VCFTOOLS
- MITHE_EXE_GNUPARALLEL
- MITHE_EXE_GATK
- MITHE_EXE_BEDOPS
- MITHE_EXE_SNPSIFT
- MITHE_EXE_SAMTOOLS

##### Workload manager environment variables
- *MITHE_NCPUS_VAR*: Must contain the name of the environment variable that indicates the number of CPU threads available for a specific job
- *MITHE_SUBMIT_CMD*: Command to submit a job. For example, "sbatch" for SLURM
- *MITHE_SUBMIT_SED*: Command that returns the ID of the job when it receives the submission command return from STDIN. For example `sed "s/Submitted batch job \(.*\)/\1/"` for SLURM
- *MITHE_SUBMIT_MUL*: Arguments to add to MITHE_SUBMIT_CMD, except the number of threads, which will be appended automatically, when submitting a multithreaded job. For example, "-N 1 -n 1 -c " for SLURM
- *MITHE_SUBMIT_PAR*: Argument to indicate the queue/partition a job should be submitted. For example, "--partition=" for SLURM 
- *MITHE_SUBMIT_DEP*: Argument to indicate the dependencies of a job when it is submitted. For example, "--dependency=afterok" for SLURM (or a different condition)
- *MITHE_SUBMIT_SEP*: Argument that separates job_ids when indicating dependencies, ":" for SLURM.
- *MITHE_MAX_TIME*: Argument to indicate the maximum time that can be allocated to a job. This will only be used for the time-consuming steps of the pipeline. We recommend setting this parameter only if the user will be using MITHE to carry out parameter optimization. For example, "-t 4-00:00" to allocate four days in SLURM.
- *MITHE_MAX_MEM*: Argument to indicate the maximum memory that can be allocated by a job. Only memory-intensive steps of the pipeline will use this amount of memory. We recommend setting this parameter only if the user will be using MITHE to carry out parameter optimization. For example, "--mem=16G" to allocate 16GB of RAM in SLURM. More than 16GB should not be needed.

##### System variables
- *MITHE_MOD_ISA*: This indicates if the local version of Environment Modules supports the is-available command or not. MITHE will work either way but will use that more advanced command if available.

## Accessory data preparation
### gnomAD
MITHE uses [gnomAD](https://gnomad.broadinstitute.org/) population allele frequency estimates. The [complete gnomad genomic VCF file is almost 500GB](https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.vcf.bgz), so we recommend to thin it keeping only the needed information for MITHE using the following commands:
```bash
gunzip -c gnomad.genomes.r2.1.sites.vcf.bgz | head -n 1000 | grep '^# ' > gnomad.genomes.r2.1.sites.onlyAFINFO.vcf
gunzip -c gnomad.genomes.r2.1.sites.vcf.bgz | perl -pe 's/\t[^\t]*;(AF=[^;]+).*$/\t$1/' >> gnomad.genomes.r2.1.sites.onlyAFINFO.vcf
```
# Usage
The main MITHE executable script is MITHE_loop.sh. This program submits a series of jobs to call SNVs and calculate the intratumor heterogeneity in a patient using trios of samples (two somatic samples, usually neoplastic, and one control sample, usually healthy tissue). MITHE calls SNVs using platypus, then filters them using a complex algorithm and additional information from public repositories, and calculates the heterogeneity between the two somatic samples. MITHE_loop.sh uses a space-separated manifest file with the following structure `Sample Normal_file SampleA_file SampleB_file [DNA]` and as many rows as patients. The filtering options are user-specified, and if several options are provided, MITHE performs this process using all possible combinations of them. The user must also specify the number of CPU threads to request multi-threaded jobs, binary flags for vcf output, list of variants output, output verbosity, and to filter out indels.

Example for parameter optimization using 8 CPU threads per parallel job, without any variant output (only the optimization summary), including indels. **WARNING**: it is not recommended to activate variant output when using MITHE for optimization due to the volume of output files (and IO operations).
```
../MITHE_loop.sh out manifest.txt params/exe_params params/filtering_params params/NAB_params params/covB_params params/PAF_params 8 0 0 0 0
```

Example for variant calling using 8 CPU threads per parallel job, with all the outputs with maximum output verbosity (all intermediate vcf files), only attending to SNVs:
```
../MITHE_loop.sh out manifest.txt params/exe_params params/filtering_params params/NAB_params params/covB_params params/PAF_params 8 1 1 2 1
```

To obtain the final list of variants after all MITHE jobs have finished, execute the MITHE_getVariants.sh, indicating the folder with MITHE results and the output file name. This can only be done when MITHE has run with output verbosity=2, and vcf output activated.

```
../MITHE_getVariants.sh out/ variants.tsv
```

## MITHE parameters
Six groups of parameters control MITHE's execution, each specified in a different input file, following the format `parameter-name,value1,value2,valuen`. Example parameter files can be found in the example directory within the repository and are explained below.

### Variant-calling parameters (exe_params)
These parameters are directly passed to the variant calling software, Platypus, in this version of the pipeline.
Example:
```
--filterReadPairsWithSmallInserts=,0
--minReads=,3
```

### Variant filtering parameters (filtering_params)
These parameters are used to filter the variants to generate the stringent variant set. Currently, the supported options are:
- -q/--qual : min quality filter
- --atoc : filter out mutations from A to C
- --atog : filter out mutations from A to G
- --atot : filter out mutations from A to T
- --ctoa : filter out mutations from C to A
- --ctog : filter out mutations from C to G
- --ctot : filter out mutations from C to T
- --gtoa : filter out mutations from G to A
- --gtoc : filter out mutations from G to C
- --gtot : filter out mutations from G to T
- --ttoa : filter out mutations from T to A
- --ttoc : filter out mutations from T to C
- --ttog : filter out mutations from T to G
- -m/--min_coverage : minimum coverage per locus
- --max_coverage: maximum coverage per locus
- -s/--min_reads_strand : minimum number of reads per strand
- -a/--min_reads_alternate : minimum number of reads for the alternative allele
- --max_reads_alternate : maximum number of reads for the alternative allele
- --min_freq_alt : min frequency of reads supporting the alternative allele

Example:
```
--qual,120
--min_coverage,0
--min_reads_strand,15
--min_reads_alternate,0
 ```

### Position filtering parameters 1: Control (NAB_params)
These parameters are used to filter somatic variants, based on information of that genomic position in the control sample. Currently, the supported filters are:

- --min_coverage : minimum coverage per locus
- --max_alternative : the maximum number of reads supporting an alternative allele
- --max_propalt : the maximum proportion of reads supporting an alternative allele

Example:
```
--min_coverage,20
--max_alternative,-1
--max_propalt,0.10
```
**NOTE**: The -1 value deactivates a filter.
**WARNING**: max_alternative and max_propralt here work differently than in covB_params, having to both be met in order to discard a variant.

### Position filtering parameters 2: B sample (covB_params)
These parameters are used to filter out variants detected as private in sample A, based on information of that genomic position in sample B. Currently, the supported filters are:

- --min_coverage : minimum coverage per locus
- --max_alternative : the maximum number of reads supporting an alternative allele
- --max_propalt : the maximum proportion of reads supporting an alternative allele

Example:
```
--min_coverage,15
--max_alternative,0
--max_propalt,0.05
```

### Population allele frequency parameters (PAF_params)
These parameters are used to filter out variants based on information on their population allele frequency. Currently, the supported filters are:

- --max_pAF: maximum population allele frequency for a variant to be kept.

Example:
```
--max_pAF,0.25
```

## Output
### Statistics and optimization
The main output of MITHE_loop.sh is the results.csv file. This comma-separated file contains many statistics at different stages of the pipeline for each case and combination of parameter values. The final similarity can be found in the filtNABcovBPAF_propU column.
### Variants
The final list of variants can be obtained using the MITHE_getVariants.sh command, as explained above.

# Citation
Fortunato A\*, Mallo D\*, Rupp SM, King LM, Hardman T, Lo J, Hall A, Marks JR, Hwang ES, Maley CC (2021) A new method to accurately identify single nucleotide variants using small FFPE breast samples. Briefings in Bioinformatics. DOI: [10.1093/bib/bbab221](https://doi.org/10.1093/bib/bbab221)
