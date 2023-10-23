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
