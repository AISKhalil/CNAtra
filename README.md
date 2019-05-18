**CNAtra** is a MATLAB-based tool for detecting large-scale segmental aneuploidies and focal amplifications/deletions of low coverage cancer cell lines. For more information please visit [CNAtra Manuscript](https://www.biorxiv.org/content/10.1101/639294v1?rss=1).

CNAtra accepts BAM/SAM files as input. It can handle both single-end and paired-end WGS reads. CNAtra pipeline has two modules: 1. Read depth calculator and 2. CNV caller. CNAtra generates many output files providing the detailed characterization of the copy number profile of the input data. It saves BED format files of both large-scale copy number variations (LCVs) and focal alterations (FAs) that can be uploaded into UCSC Genome Browser. In addition, CNAtra incorporates a visual platform that allows manual review of the CNV profile and its accessory information such as the mappability score. The detailed description of CNAtra inputs, parameters, methods, and outputs are provided in **CNAtra_User_Guide.pdf"** file.

**CNAtra output at a glance**

User can run CNAtra tool and get the CNV profile quickly using the default parameters if the BAM file is available. We provided the simulated data where we artificially incorporated LCVs and FAs in Chr3 of CHP-212 cell line data as an example for testing CNAtra installation. Input BAM is available as "CHP212_chr3_Artificial.bam" file under "CNAtraInput" folder. CNAtra outputs were provided in "CNAtraResults" folder (default output directory). After cloning CNAtra git repository, user can run the following commands using MATLAB to get the CNV profile.

​			*>> CNAtraDirectory = './CNAtraTool';*

​			*>> addpath(CNAtraDirectory);*

​			*>> inputFile = './CNAtraInput/CHP212_chr3_Artificial.bam';*

​			*>> CNAtraObj = CNAtra(inputFile, CNAtraDirectory);*

​			*>> CNAtraObj.RDcalculator;*

​			*>> CNAtraObj.CNVcaller;*
