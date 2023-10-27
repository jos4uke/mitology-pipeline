

# mitology-pipeline

Pipeline for plant organites (chloroplast and mitochondrion) assembly

## Badges

[![GPLv3 License](https://img.shields.io/badge/License-GPL%20v3-yellow.svg)](https://opensource.org/licenses/)

## Dependencies  

The following dependencies need to be installed manually:
* [bash-common](https://github.com/jos4uke/bash-common.git)  
* [log4sh](https://github.com/kward/log4sh)  

The following dependencies can be installed using [conda](https://docs.conda.io/projects/miniconda/en/latest/) (recommended):
* [bgzip](http://www.htslib.org/)    
* [bwa](https://github.com/lh3/bwa)    
* [bcftools](http://www.htslib.org/)  
* [khmer](https://github.com/dib-lab/khmer)  
* [MetaVelvet](http://metavelvet.dna.bio.keio.ac.jp/)
* [Velvet](https://github.com/dzerbino/velvet)      
* [nucmer](https://github.com/mummer4/mummer)  
* [quast](https://github.com/ablab/quast)  
* [samtools](http://www.htslib.org/)  
* [tabix](http://www.htslib.org/)  

### install bash-common

Clone the `bash-common` repository then install the bash library (tag: mitology) system-wide into `/usr/local` under `share` directory (need to have at least sudoer privileges):    
  
```bash
git clone https://github.com/jos4uke/bash-common.git 
cd bash-common
sudo bash install.sh mitology /usr/local

```

### install log4sh

Clone the `log4sh` repository into `/usr/loca/share` (need at least to have sudoer privileges):  

```bash
cd /usr/local/share
sudo git clone https://github.com/kward/log4sh.git

```

## Installation

```bash
cd /usr/local/share
sudo git clone https://github.com/jos4uke/mitology-pipeline.git
cd mitology-pipeline

```

### File content

```
# tree
├── bin     # contains the main pipeline script
│   └── mitology-pipeline.sh
├── COPYING
├── README.md
└── share   # contains the main pipeline script
    └── mitology-pipeline
        ├── etc     # contains the pipeline configuration file to be copied by the end user
        │   └── mitology-pipeline_user.config
        ├── lib     # contains specific bash library functions
        │   ├── mitology-alignments_lib.inc
        │   └── mitology-pipeline_lib.inc
        └── scripts     # helper scripts
            ├── plot_hashcount.sh
            ├── R
            │   ├── length-weigthed_kmer_coverage_hist.rplot.R
            │   └── plot_hashcount_hist.R
            ├── run_meta-velvetg.sh
            └── split_files.sh

7 directories, 12 files

```

## Usage guide  

**Steps:**
1. get a copy of the configuration file  
2. update the configuration file  
3. run the pipeline  


### get a copy of the configuration file

```bash
cp share/mitology-pipeline/etc/mitology-pipeline_user.config .

```

### update the configuration file  

The configuration file contains different sections:  
* paths  
* genome_alias  
* sample  
* khmer_load_into_counting  
* khmer_filter_abund  
* contig_assembler  
* scaffolder  
* velveth  
* velvetg  
* meta_velvetg   
* quast  
* bwa_aln  
* bwa_sampe  
* filtering  
* samtools_view  
* samtools_mpileup  
* bcftools_view  
* bgzip  

In the `[paths]` section, you need to update path to all resources listed as `<key>=<value>` pair.  
- `GENOMES_BASE_PATH` and `INDEXES_BASE_PATH` are the root directory path to genomes files and tool index respectively.  
- `BWA_INDEXES` and `SAMTOOLS_INDEXES` are the location of the bwa and samtools indexes directories respectively.  
- `bcftools`, `bgzip`, `bwa`, `khmer_*`, etc. are the paths to the corresponding tools. They should be reachable from your `PATH` else please indicate the absolute path. It is recommended to install all these tools using [conda](https://docs.conda.io/projects/miniconda/en/latest/).   
  
In the `[genome_alias]`, you need to update all prefixes/aliases to reference genomes used notably for indexes listed as `<key>=<value>` pair. `ref` and `ref_short` aliases refer to the host genome, whereas `ref_mito*` and `ref_chloro*` aliases refer to the mitochondrial and chloroplast genomes or annotations (gff) respectively. 

In the `[sample]` section, you need to update the sample alias and paths listed as `<key>=<value>` pair. `name_alias` will be used to identify sample in the pipeline. And `seqfile_parent_dir` is the parent directory location of the paired-end reads files specified by `seqfile_{R1,R2}` keys.  
  
The other sections refer to options specific to pipeline steps and tools, listed again as `<key>=<value>` pair. Please refer to the tool documentation to update the options values.  
  
### run the pipeline

Run the pipeline using the updated config file, and output results into `results` directory.  

```bash
bash bin/mitology-pipeline.sh -c mitology-pipeline_user.config -o results 

``` 

Please checkout the usage help with  `bash mitology-pipeline.sh --help`

#### ouput files

todo: 
* describe here all output files
* need some test dataset 

## Authors

- [@jos4uke](https://www.github.com/jos4uke)


## Acknowledgements

IJPB Bioinformatics team

## License

[![GPL-3.0-or-later](assets/img/gplv3-or-later.png)](ttps://www.gnu.org/licenses/gpl-3.0.fr.html)  
see `COPYING`


