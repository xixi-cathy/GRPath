<!-- ABOUT GRPath -->

## About GRPath

**GRPath** is a computational framework that can calculate the genotype-to-phenotype gene regulation paths in the form of **"causal genomic region-causal variant-causal eGene-noteworthy cell type-disease phenotype"** for diseases.

<img src=".\figure\fig.jpg" alt="overview" style="zoom:10%;" />



<!-- AUTHOR -->

## Author

Xi Xi xix19@mails.tsinghua.edu.cn


  
<!-- GETTING STARTED -->

## Getting started

### Software prerequisites

#### Environment

```
python >= 3.6.10
R >= 3.6.1
```



#### Install Python packages

```shell
pip install math
pip install random
pip install pandas
pip install numpy
pip install sklearn
pip install scipy
pip install multiprocessing
pip install functools
pip install collections
```



#### Install R packages

```R
install.packages("dplyr")
install.packages("tidyr")
```
  

### Data preparation

You should first decide the disease you would like to study, then prepare the following data:

- **Whole-genome sequencing (WGS)**,  **RNA-seq** data of disease-relevant tissues and **disease phenotype** (disease/control) from the <font color='red'>**same**</font> individual.

  If accessible, protected data from the [GTEx](https://gtexportal.org/home/protectedDataAccess) project are more than appropriate.

  Here, if the disease phenotype information are not that detailed, you may relax the disease to a broader class that covers the one you would like to study. For example, we would like to study heart failure caused by coronary heart disease (cHF), so we used "cardiovascular disease" to label each individual.

- **eQTLs of the disease-relevant tissues**.

  You may obtain tissue-eQTLs from the [GTEx data download](https://gtexportal.org/home/datasets) website.

- **GWAS summary statistics** of the disease you would like to study.

  The disease must be the same as the one you used to label each individual. For example, cardiovascular disease.

- **scRNA-seq** data of disease and control samples.

  The scRNA-seq data are not necessarily from the same samples that provide WGS, RNA-seq and disease phenotype, but the disease cells should be of the exact disease you would like to study. For example, cHF.  



<!-- USAGE -->
## Usage

### 1. Calculate the chromatin openness scores

First, you may refer to the method proposed by [Li et al.](https://www.pnas.org/content/117/35/21364) to predict **(disease-relevant-) tissue-specific** openness scores of pre-defined regulatory elements (REs) corresponding to variants in the donor population. This step should produce 3 files:

- **variant_openness.txt**: a variant-by-donor matrix that denotes the openness scores of each variant predicted from personal genomes.

  The first 8 columns should be the RE and variant information, which are RE chromosome, RE starting site, RE ending site, variant chromosome, variant starting site, variant ending site, variant reference allele and variant alternate allele.

- **variant_allele.txt**: a binary variant-by-donor matrix (the same variants as in "variant_openness.txt") that records the genetic variations of each sample (0: reference allele / 1: alternate allele).

  The first 8 columns should be the same as in "variant_openness.txt".

- **RE_openness.txt**: an RE-by-sample matrix that denotes the openness scores predicted from reference genomes of all REs that contain the variants in "variant_openness.txt". Note: a donor may provide multiple samples.

  The first 3 columns should be the RE information, which are RE chromosome, RE starting site and RE ending site.  



### 2. Predict causal genomic regions, causal variants and corresponding causal eGenes

Next, you will get some preliminary results, which are part of the regulation path in the form of "causal genomic region-causal variant-causal eGene" of the disease you are interested in. Input files should include:

- **variant_openness.txt**, **variant_allele.txt** and **RE_openness.txt** from the previous step.
- **GWAS.txt**: location of the GWAS SNPs, which include variant name, chromosome name, and chromosome position.
- **phenotype.txt**: donor disease phenotype information. The 4 columns must be donor name, the number of sample in "RE_openness.txt", the number of sample in "variant_openness.txt", and disease label (0: control / 1: disease).
- **variant_pos.txt**: location of variants in the form of "Chromosome_Position_ReferenceAllele_AlternateAllele_b37". For example, "chr2_136422784_C_T_b37" for "rs3863014".
- **eqtls_200kb.txt**: tissue-eQTLs within the 200 kb region around each GWAS SNP. The 3 columns should be region central SNP (GWAS SNP), variant name and corresponding eGene.



After preparing all these data, run:

```sh
bash predict_region_variant_egene.sh
```



The output files are:

- **alter_pathogenic_region_variant_egene.csv**: alternate-allele-pathogenic regions and causal variants within and corresponding causal eGenes.
- **ref_pathogenic_region_variant_egene.csv**: reference-allele-pathogenic regions and causal variants within and corresponding causal eGenes.
- A "processing" folder that contains some intermediate products.

Checking "alter_pathogenic_region_variant_egene.csv" and "ref_pathogenic_region_variant_egene.csv", you may find some interesting biological stories.  



### 3. Predict noteworthy cell types and get the complete gene regulation paths

Now, if you want to dig further, you could incorporate scRNA-seq data of disease and control samples, and extend the regulation paths to cell type level. You should prepare the following files:

- **single_cell_data.txt**: a gene-by-cell matrix. Values in the matrix are **normalized gene expression levels**. So, you should process your raw scRNA-seq data beforehand.
- **single_cell_label.txt**: disease state (0: control / 1: disease) and cell type of each cell in the gene-by-cell matrix.



After preparing these files, run:

```sh
bash predict_regulation_path.sh
```



The output files are:

- **regulation_path.csv**: final regulation paths in the form of "causal genomic region-causal variant-causal eGene-noteworthy cell type-disease phenotype".
- cell_type_auc_mean.csv and cell_type_auc_std.csv in the "processing" folder as intermediate products.  



<!-- EXAMPLE -->

## Example

You may download the demo_data and run:

```sh
bash predict_region_variant_egene.sh
bash predict_regulation_path.sh
```

to find the regulation paths for heart failure caused by coronary heart disease (cHF).  



<!-- ACKNOWLEDGEMENTS -->

## Acknowledgements

* Li, W., Duren, Z., Jiang, R., Wong, W.H., 2020. A method for scoring the cell type-specific impacts of noncoding variants in personal genomes. Proc Natl Acad Sci USA 117, 21364–21372. [(Link)](https://doi.org/10.1073/pnas.1922703117)
* Wang, L., Yu, P., Zhou, B., Song, J., Li, Z., Zhang, M., Guo, G., Wang, Y., Chen, X., Han, L., Hu, S., 2020. Single-cell reconstruction of the adult human heart during heart failure and recovery reveals the cellular landscape underlying cardiac function. Nat Cell Biol 22, 108–119. [(Link)](https://doi.org/10.1038/s41556-019-0446-7) 
* Skinnider, M.A., Squair, J.W., Kathe, C., Anderson, M.A., Gautier, M., Matson, K.J.E., Milano, M., Hutson, T.H., Barraud, Q., Phillips, A.A., Foster, L.J., La Manno, G., Levine, A.J., Courtine, G., 2021. Cell type prioritization in single-cell data. Nat Biotechnol 39, 30–34. [(Link)](https://doi.org/10.1038/s41587-020-0605-1)

