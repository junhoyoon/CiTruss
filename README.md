# CiTruss
This software implments CiTruss, a statistical framework called CiTruss for simultaneously learning a gene network and *cis*-acting and *trans*-acting eQTLs that perturb this network, given population allele-specific expression and SNP data.

# Installation

## Linux

Go to ./CiTruss folder and run the following command to install the software.

```bash
./make
```

## Mac OSX

Depending on the MacOSX version and default compiler, the location of some of the library #include statements in Mega-sCGGM, a software required to run CiTruss. Just always using the llvm/clang compiler installed via homebrew puts the libraries in the same place. For this, you will first need to install [Homebrew](https://brew.sh/), and then use it to install llvm:

```bash
brew install llvm
```

You will also need to make sure you have the latest versions of XCode command-line tools installed.

Then, go to ./CiTruss folder and run
```bash
./make
```

# Usage & Demo

## Sum model

Usage
```bash
./citruss_sum num_samples num_genes num_snps total_gene_expression genotype_allele1 genotype_allele2 regularization_paramemter_for_network regularization_paramemter_for_eQTL -o output_path
```

Demo
```bash
./citruss_sum 200 40 100 ./demo/Ysum.txt ./demo/Xm.txt ./demo/Xp.txt 0.1 0.1 -o ./
```

## Difference model

Usage
```bash
./citruss_diff num_samples num_genes num_snps gene_expression_allele1 gene_expression_allele2 genotype_allele1 genotype_allele2 gene_info snp_info regularization_paramemter_for_cis_eqtl -o output_path -r 1 --init-Psi eqtl_file
```

Demo
```bash
./citruss_diff 200 40 100 ./demo/Ym.txt ./demo/Yp.txt ./demo/Xm.txt ./demo/Xp.txt ./demo/gene_info.txt ./demo/snp_info.txt 0.1 -o ./demo_o/ -r 1 --init-Psi ./demo/eqtl_init.txt
```
