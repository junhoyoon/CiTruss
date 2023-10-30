# CiTruss
This software implments CiTruss, a statistical framework for simultaneously learning a gene network and *cis*-acting and *trans*-acting eQTLs that perturb this network, given population allele-specific expression and SNP data. CiTruss consists of 
- [sum model](#sum-model) that takes total gene expression and phased genotype data to learn a gene network and eQTLs together with target genes
- [difference model](#difference-model) that takes allele-specific gene expression, phased genotype data, and eQTLs from the sum model to learn *cis*-acting eQTLs together with target genes.

Below, we provide instructions on
- [how to install CiTruss](#installation)
- [format gene expression and SNP data](#data-formatting)
- [run CiTruss](#usage--demo)
- [understand outputs from CiTruss](#understanding-the-citruss-outputs).

CiTruss was introduced in the following preprint:
> Jun Ho Yoon and Seyoung Kim. **Learning gene networks under SNP perturbation using SNP and allele-specific expression data**. bioRxiv 2023.10.23.563661; doi: https://doi.org/10.1101/2023.10.23.563661

## System Requirements

CiTruss was compiled and tested on the following environments:

- Linux: CentOS 7.4 with g++ (GCC) 4.8.5
- MacOS:

CiTruss also uses [Metis 5.1.0](http://glaros.dtc.umn.edu/gkhome/metis/metis/download), which is already included in [./metis-5.1.0](/metis-5.1.0).

## Installation

### Linux

Inside ./CiTruss folder, run the following command to install CiTruss. It will take less than a minute to compile CiTruss.

```bash
./make
```

This will generate `citruss_sum` for [CiTruss sum model](#sum-model) and `citruss_diff` for [CiTruss difference model](#difference-model).

### Mac OSX

Depending on the MacOSX version and default compiler, the location of some of the library #include statements in Mega-sCGGM, a software required to run CiTruss. Just always using the llvm/clang compiler installed via homebrew puts the libraries in the same place. For this, you will first need to install [Homebrew](https://brew.sh/), and then use it to install llvm:

```bash
brew install llvm
```

You will also need to make sure you have the latest versions of XCode command-line tools installed. Then, go to ./CiTruss folder and run the following command to install CiTruss. It will take less than a minute to compile CiTruss.

```bash
./make
```

This will generate `citruss_sum` for [CiTruss sum model](#sum-model) and `citruss_diff` for [CiTruss difference model](#difference-model).

## Data Formatting

CiTruss takes gene expression data (total and allele-specific expression) and phased genotype data. The difference model additionally requires eQTL, gene, and SNP information. Below are instructions on how to format these files.

### Gene expression

Gene expression data should be a text file that includes a space-separated matrix of normalized gene expression levels in *(number of genes, number of samples)* shape. The sum model takes total gene expression (an example in [/demo/Ysum.txt](./demo/Ysum.txt)), and the difference model takes allele-specific expression (examples in [/demo/Y1.txt](./demo/Y1.txt) and [/demo/Y2.txt](./demo/Y2.txt)). When allele-specific expression is not available for certain genes or samples, the entries should be "NaN". For normalization, we used *log2* of TPM (Transcripts Per Million).


### Genotype

Genotype data should be a text file that includes a space-separated matrix of phased genotype information in *(number of SNPs, number of samples)* shape. Both the sum and difference models take phased genotype data (examples in [/demo/X1.txt](./demo/X1.txt) and [/demo/X2.txt](./demo/X2.txt)).


### eQTL information

eQTL information should be a text file that includes indices of SNPs along with the indicies of their target genes in each line. The indicies of SNPs and genes should match the row indicies of SNPs and genes in allele-specific expression and phased genotype files. A score metric for each pair of eQTL and target gene should also be included. Lastly, the first line should contain the number of SNPs, genes, and eQTLs, each separated with a space (an example in [/demo/eqtl_init.txt](./demo/eqtl_init.txt)).

For score metrics, an output from CiTruss sum model can be directly used. If eQTLs were obtained with other softwares that outputs p-values, for instance, a negative log of p-value can be used.

For example, the first two lines of eQTL information for 100 SNPs, 40 genes, and 200 eQTLs should look like
```
100 40 200
1 15 -0.0746932222127
```

### Gene & SNP information

Gene information should be a text file that includes the chromosome, gene start position, and gene end position for each gene in each line. The row indicies should match those in allele-specific expression files (an example in [/demo/gene_info.txt](./demo/gene_info.txt)). Similarly, SNP information should be a text file that includes the chromosome and position for each SNP in each line (an example in [/demo/snp_info.txt](./demo/snp_info.txt)).


## Usage & Demo

CiTruss consists of a **sum model** that takes total gene expression and phased genotype data to learn a gene network and eQTLs together with target genes and a **difference model** that takes allele-specific gene expression, phased genotype data, and eQTLs from the sum model to learn *cis*-acting eQTLs together with target genes. Below, we provide usage and demo run for each of the sum and difference models.

### Sum model

```
Usage:
./citruss_sum num_samples num_genes num_SNPs Ysum_file X1_file X2_file reg_V reg_F [options]\n

required arguments:
    num_samples: number of samples
    num_genes: number of genes
    num_SNPs: number of SNPs
    Ysum_file: file path for total gene expression
    X1_file: file path for phased genotype of haplotype 1
    X2_file: file path for phased genotype of haplotype 2
    reg_V: regularization parameter for V [gene network]
    reg_F: regularization parameter for F [eQTLs]

options [-flag <option> (default option)]:
    -o <output_prefix> (./): prefix for output files
    --init-V <V0_filename> (none): file path for initial V [gene network]
    --init-F <F0_filename> (none): file path for initial F [eQTLs]
    -v <verbose> (1): show optimization information or not (1 or 0)
    -i <max_iters> (10): max number of outer iterations
    -s <sigma> (1e-4): backtracking termination criterion
    -q <tol> (1e-2): tolerance for terminating outer loop
    -j <obj_tol> (1e-13): CG tolerance for calculating objective function
    -g <grad_tol> (1e-10): CG tolerance for calculating gradient
    -h <hess_tol> (1e-8): CG tolerance for calculating hessian
    -m <memory_usage> (32000): memory capacity in MB
    -n <threads> (16) : set the max number of threads
    -r <refit> (false): update gene network and eQTLs without adding new entries from init files
```

A demo run for sum model:
```bash
./citruss_sum 200 40 100 ./demo/Ysum.txt ./demo/X1.txt ./demo/X2.txt 0.1 0.1 -o ./demo/
```

The demo run should finish in a few seconds. The sum model generates two outputs, `V.txt` for gene network and `F.txt` for eQTLs. Expected outputs from the demo run can be found in [/demo](./demo).


### Difference model


```
Usage:
./citruss_diff num_samples num_genes num_SNPs Y1_file Y2_file X1_file X2_file gene_info snp_info reg_Psi [options]

required arguments:
    num_samples: number of samples
    num_genes: number of genes
    num_SNPs: number of SNPs
    Y1_file: file path for allele-specific gene expression of haplotype 1
    Y2_file: file path for allele-specific gene expression of haplotype 2
    X1_file: file path for phased genotype of haplotype 1
    X2_file: file path for phased genotype of haplotype 2
    gene_info: file path for gene information
    snp_info: file path for SNP information
    reg_Psi: regularization parameter for Psi [cis-acting eQTLs]

options [-flag <option> (default option)]:
    -o <output_prefix> (./): prefix for output files
    --init-Gamma <Gamma0_filename> (none): file path for initial Gamma [variance for allele-specific expression of each gene]
    --init-Psi <Psi0_filename> (none): file path for initial Psi [cis-acting eQTLs]
    -v <verbose> (1): show optimization information or not (1 or 0)
    -i <max_iters> (10): max number of outer iterations
    -s <sigma> (1e-4): backtracking termination criterion
    -q <tol> (1e-2): tolerance for terminating outer loop
    -j <obj_tol> (1e-13): CG tolerance for calculating objective function
    -g <grad_tol> (1e-10): CG tolerance for calculating gradient
    -h <hess_tol> (1e-8): CG tolerance for calculating hessian
    -m <memory_usage> (32000): memory capacity in MB
    -n <threads> (16) : set the max number of threads
    -r <refit> (false): update cis-acting eQTLs without adding new entries from init file
```


A demo run for difference model:
```bash
./citruss_diff 200 40 100 ./demo/Y1.txt ./demo/Y2.txt ./demo/X1.txt ./demo/X2.txt ./demo/gene_info.txt ./demo/snp_info.txt 0.1 -o ./demo/ -r 1 --init-Psi ./demo/eqtl_init.txt
```

The demo run should finish in a few seconds. The difference model generates two outputs, `Psi.txt` for *cis*-acting eQTLs and `Gamma.txt` for variance of allele-specific expression. Expected outputs from the demo run can be found in [/demo](./demo).


## Understanding the CiTruss Outputs


