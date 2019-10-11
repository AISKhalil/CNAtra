# CNAtra: An analytical and visualization tool for hierarchical CNV discovery of large-scale and focal copy number alterations in low-coverage cancer genomes

`CNAtra` is a MATLAB-based tool that accepts BAM/SAM files as input. It can handle both single-end and paired-end WGS reads. `CNAtra` pipeline has two modules: 1. Read depth calculator and 2. CNV caller. `CNAtra` generates many output files providing the detailed characterization of the copy number profile of the input data. It saves BED format files of both large-scale copy number variations (LCVs) and focal alterations (FAs) that can be uploaded into UCSC Genome Browser. In addition, `CNAtra` incorporates a visual platform that allows manual review of the CNV profile and its accessory information such as the mappability score. The detailed description of `CNAtra` inputs, parameters, methods, and outputs are provided in **"CNAtra_User_Guide.pdf"** file.

**. For a full description of the method and applications, please visit [CNAtra Manuscript](https://www.biorxiv.org/content/10.1101/639294v1).**

## CNAtra output at a glance

User can run `CNAtra` tool and get the CNV profile quickly using the default parameters if the BAM file is available. We provided the simulated data where we artificially incorporated LCVs and FAs in Chr3 of CHP-212 cell line data as an example for testing `CNAtra` installation. Input BAM is available as "CHP212_chr3_Artificial.bam" file under "CNAtraInput" folder. `CNAtra` outputs were provided in "CNAtraResults" folder (default output directory). 

### Download CNAtra from github

```bash
cd ~
git clone https://github.com/AISKhalil/CNAtra.git
```

### Use MATLAB to get the CNV profile of BAM file

```bash
>> CNAtraDirectory = './CNAtraTool';
>> addpath(CNAtraDirectory);
>> inputFile = './CNAtraInput/CHP212_chr3_Artificial.bam';
>> CNAtraObj = CNAtra(inputFile, CNAtraDirectory);
>> CNAtraObj.RDcalculator;
>> CNAtraObj.CNVcaller;
```
