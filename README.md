# xQTLImp
## Introduction
xQTLImp is an open source software that implements imputation of xQTL(such as eQTL, mQTL, haQTL et. al) statistics across the genome. xQTLImp accepts xQTL summary statistics (i.e., Z statistics) without the need of individual-level genotypes and molecular traits (such as gene expression profiles), and could accurately impute novel xQTL associations from known associations based on linkage disequilibrium (LD) among variants. Specifically, it models the statistics of variants associated with the same molecular trait by using a multivariate Gaussian model, and infers missing association statistics from nearby known association statistics by leveraging LD among variants (See Methods in the paper). Using multiple real datasets, we demonstrated that 1) xQTLImp can impute missing xQTL statistics with high accuracy, and 2) xQTLImp can effiectively reduce the lower bound of MAF in xQTL studies, leading to the discovery of novel xQTL signals to further enhance xQTL discoveries.

xQTLImp requires as compulsory input data: (1) genotype reference panel (such as 1000G or HapMap) for calculating LD correlations; (2) xQTL summary statistics (Z statistics) measured from associations between genotypes and molecular traits, where variants are identified by genomic coordinates and ref/alt alleles; (3) molecular trait annotation. The output is a list of variant-molecular trait associations, with a imputation flag (1 for novel) and imputation quality score, r2pred from (0,1). 

xQTLImp could handle various kinds of molecular traits (not limited to gene expression, DNA methylation or histone acetylation) that can be physically mapped onto the genome. Genomic variants including SNV and small Indel (Insertion/deletion) are both supported. xQTLImp allows users to specify chromosome and MAF of variants, discard certain regions (such as HLA region) and run in multiple threads. Regarding the performance, xQTLImp was designed in speed and memory efficient way, and it took 3~4 hours and <4GB memory for the genome-wide imputation on a single-cell eQTL datasets using 20 threads (see Suppl. Notes in the paper), which could be easily distributed on PC and server.

The source code and sample data are freely available for download at current webpage. The software is released under GNU GPL license (version 3). xQTLImp is implemented in C++ and has been successfully tesed in Linux and Mac OS X platforms.
</br>

##  Building xQTLImp
 To read gzipped VCF files, **zlib** package is required, and to install:
```bash
sudo apt-get install zlib1g-dev
```
Then compiling:
```bash
cd ./src/
make #Any C++11 compiler should work.
```
## Usage
### Requirements for input files.
The following files and format are required as input:
#### 1. Genome reference panel in gzipped VCF format. 
For example, the 1000G human genome reference panel [(Available here)](http://www.internationalgenome.org/data/), or HapMap3 reference panel [(Available here)](https://www.sanger.ac.uk/resources/downloads/human/hapmap3.html).</br>

```bash
# Codes for downloading 1000G Chr1 to Chr22 VCF files.
prefix="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr" ;
suffix=".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" ;
for chr in {1..22} ; do
    wget $prefix$chr$suffix  $prefix$chr$suffix.tbi ;
done
```
The genotype reference VCF files should be separated by chromosomes, and named in the format of 'chr.*N*.XXX.vcf.gz' within a folder. *N* represents for chromosome number (interger), for example 1~26 (X->23, Y->24, XY (Pseudo-autosomal region of X) ->25, MT (Mitochondrial) ->26); XXX represents for any user defined string. Other domain should be freezed in required format. </br>
</br>

*Notes:* To get better performance and accuracy, we recommend the users to preprocess the genotype reference VCF files. For example, keeping only samples of EUR population if the input xQTL statistics are derived from subjects of European ancestry, or/and filtering rare and non-biallelic variants since most of the current xQTL studies still focus on common variants due to limited sample size. The example source codes with *[VCFtools](https://github.com/vcftools/vcftools)* are as follows:

```bash
for chr in {1..22} ; do
vcftools --gzvcf $prefix$chr$suffix \
--maf 0.01                          \ # keep common variants with MAF > 0.01
--min-alleles 2 --max-alleles 2     \ # keep biallelic variants
--keep $EUR_file                    \ # keep EUR population samples, not provided.
--remove-filtered-all --recode --stdout | gzip -c > $out_prefix$chr$out_suffix;
done
```

#### 2. Molecule annotation file
Molecule annotation file gives the physical position of each molecule on reference genome. This annotation file should start with a column name line that contains at least three columns – molecular_ID, start_pos, end_pos ...(optional columns), followed by lines of data entries. Each field of data entries must be separated by tab or white spaces. Header names can be different with the example.</br>
##### Example:

`molecular_ID`| `start_pos`| `end_pos`
------------------|-------|-------
ENSG00000223972.4	| 11869 | 14412
ENSG00000227232.4	| 14363	| 29806
ENSG00000243485.2	| 29554	| 31109
ENSG00000237613.2	| 34554	| 36081
......</br>


#### 3. xQTL summary statistics
This file includes summary statistics associated with pairs of variants and molecular traits. The xQTL file should start with a line that contains at least 6 columns (column names can be different) – chromosome , molecular_ID, variant_pos , ref_allele , alt_allele , z_statistics ...(optional columns), followed by lines of data entries . Each field of data entries must be separated by tab or white spaces. Data entries should be **sorted** at least by chromosome number in increasing order, and records with same molecular_ID should be grouped together. We recommend users to sort the xQTL file by chr, molecular_start_pos, molecular_ID, and variant_pos in increasing order prior imputation. Users can also specify ```--sort=TRUE``` if their input xQTL file is not sorted properly.</br>
##### Example:
`chromosome` | `molecular_ID` | `variant_pos` | `ref_allele` | `alt_allele` | `z_statistics`
--|--|--|--|--|--
1 | ENSG00000223972 | 13417 | G | C    | 1.5
1 | ENSG00000227232 | 17559 | A | AGCC | 2.6
1 | ENSG00000227232 | 54421 | G | A    | -1.0
......</br>

*Notes:* </br>
* xQTLImp only handle cis-xQTL associations because of the nature of LD, so variant and molecule should be on the same chromosome. </br>
* The coordinates of variants and molecules should be always in same genome build (e.g. hg19).
* The Z statistic should be computed as the effect of the same type of allele (Ref or Alt) for all data entries (i.e. The effect allele is Ref allele or Alt allele).  And the effect allele of imputed variants will be the same type of effect allele as the inputs.
* The sign of Z statistics will be converted by xQTLImp if the user provided <ref, alt> alleles are in opposite to that in the reference panel.
* Data entries (input lines) would not be used for imputation if xQTLImp cannot locate the variant on reference panel, but they will still be written into outputs without loss of information. 


### xQTLImp parameters：
The parameters can be specified by short options (e.g. -h) or long options (e.g. --help), which have same effects. Please see the parameters and explanations as follows:
```bash
xQTLImp
-h, --help           null           # null, to display this usage.
-x, --xQTL           file_path      # string, the file path of xQTL summary statistics.
-m, --molecule       file_path      # string, the file path of molecule annotation file.
-v, --VCF            folder_path    # string, the folder path of genome reference panel, such as 1000G VCF files.
-o, --output         folder_path    # string, the folder path of output results. 
-c, --chr            chromosome     # int, specify which chromosome will be imputed.
-t, --num_threads    num_threads    # int, number of threads, 1 in default.
-s, --sort           TRUE/FALSE     # boolean, sort the xQTL summary statistics by chromosome, molecular_ID, and variant_pos in increasing order prior imputation (required), FALSE in default.
-e, --exclude        chr:start-end  # int:int-int, specify a genome region in which variants will be ignored during imputation process.
-b, --exclude_file   file_path      # string, multiple genome regions user want to mask during imputation process.
-f, --MAF_cutoff     MAF_cutoff     # double, minimum MAF threshold for variants in genome reference panel, 0.01 in default.
-l, --lambda_value   lambda_value   # double, a constant value used to added with var-covariance matrix to gurantee the matrix is invertible, 0.1 in default. 
-w, --window_size    window_size    # int, Window size N, +-N/2 apart from molecular center pos, in base pair, 500000bp in default.
```
*Notes:* 
* ```-e or --exclude``` is useful for users to exclude genome regions they want to ignore, such as the complex human HLA region (6:25000000-35000000 in hg19), which will tremendously slow down the imputation process because of high density of genomic markers. 
* Combining ```-c or --chr``` with ```-t or --num_threads``` could greatly save time for users by taking advantage of HPC. For example, on a slurm based cluster, the following script named *Cluster_sample.sh* can be submitted onto multiple nodes with multiple threads:

```bash
cat Cluster_sample.sh  # in shell
#!/bin/bash
#SBATCH -c 16                   # Number of cores
#SBATCH -t 2:00:00              # Runtime
#SBATCH -p short                # Partition (queue) to submit to
#SBATCH --mem-per-cpu=1G        # memory needed (memory PER CORE)
chr=$1
Path_to_xQTLImp/src/xQTLImp             \
-x your_path/Input_xQTL_file            \
-m your_path/molecular_annotation_file  \
-v your_path/1000G_ref_panel/           \
-o your_output_folder_path/             \
-e 6:25000000-35000000                  \
-c $chr                                 \
-t 16

# shell code
for chr in {1..22}}; do
	sbatch Cluster_sample.sh $chr
done
```

### Running sample data
The sample of summary statistics are generated from eQTL summary data of CD4+ T cells from a single-cell eQTL study by Monique et. al, published in [Nature Genetics](https://www.nature.com/articles/s41588-018-0089-9) in 2018. The input eQTL summary file includes 24,846 variant-gene associations representing 100 genes and 5106 variants on chromosome 1. The sample of reference panel is generated from 1000G phase3 (EUR population) with variants having MAF>0.01. 198,337 associaitons including the inputs will be generated by xQTLImp. The whole process will take approximately 5~10 mins.


```bash
cd ./src ; make; cd .. # compliling under src folder
mkdir sample_output    # creat a new folder for output
# Running xQTLImp
./src/xQTLImp                                  \
-m ./sample/sample_gene_annotation.txt         \
-x ./sample/sample_eQTL_summary.txt            \
-v ./sample/                                   \
-o ./sample_output/                            \
-t 1                                           \
-w 500000
```
### xQTLImp output file format
*N* imputation results corresponding with *N* chromosomes will be created in the output folder.</br> 
#### Example:
`Chr` | `Molecular_ID` | `Molecular_Start` | `Molecular_End` | `Variant_ID` | `Variant_pos` | `Variant_Ref` | `Variant_Alt` | `Z-Statistic` | `R2pred` | `Imputation_flag`
--|--|--|--|--|--|--|--|--|--|--
1 | ENSG00000231709 | 521369 | 523833 | 1_13417_C_CGAGA | 13417 | C | CGAGA | -0.338210 | 0.98 | 1
1 | ENSG00000231709 | 521369 | 523833 | 1_17559_G_C | 17559 | G | C | 1.605512 | 1.000000 | 0
1 | ENSG00000231709 | 521369 | 523833 | 1_54421_A_G | 54421 | A | G | 1.069012 | 1.000000 | 0
......

## Reference
During the implementation of xQTLImp, we referenced the following works:
* Pasaniuc, B. et al., 2014. Fast and accurate imputation of summary statistics enhances evidence of functional enrichment. Bioinformatics, 30(20), pp.2906–2914.
* Kwan, J.S.H. et al., 2016. FAPI: Fast and Accurate P-value Imputation for genome-wide association study. European Journal of Human Genetics, 24(5), p.761.
* Han, B., Kang, H.M. & Eskin, E., 2009. Rapid and accurate multiple testing correction and power estimation for millions of correlated markers. PLoS genetics, 5(4), p.e1000456.
* van der Wijst, Monique GP, et al. "Single-cell RNA sequencing identifies celltype-specific cis-eQTLs and co-expression QTLs." Nature genetics 50.4 (2018): 493.


## Bug shooting

* Mac user might occur the following error during compiling:
```
xcrun: error: invalid active developer path (/Library/Developer/CommandLineTools), missing xcrun at: /Library/Developer/CommandLineTools/usr/bin/xcrun
``` 
which means you need to install proper C++ compiling environment. Installing XCode can fix this problem, open a Terminal and run this command:
```
xcode-select --install
```

* Segmentation fault:

This error might be due to various reasons. Please first double check if the parameters you specified satisfy our requirements. For example, -v parameter requires a folder path rather than a file path. And VCF files should be named in the format as "chr.1.user_defined_string.vcf.gz". This is because when multiple chromosomes exist in your xQTL summaries, which is the most common case, xQTLImp will automatically search the folder for VCF files under the naming scheme.
If you have confidence in your inputs, please send us an email with your commands and light sample data.

