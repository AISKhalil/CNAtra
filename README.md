**CNAtra** is a MATLAB-based tool for detecting large structural-variations and focal amplifications/deletions of low coverage cancer cell line. For more information please visit [CNAtra Manuscript](https://www.biorxiv.org/).

CNAtra accepts BAM/SAM files as input. It can handle both single-end and paired-end WGS reads. CNAtra pipeline has two main modules: 1. Read depth calculator and 2. CNV caller. CNAtra generates many output files providing the detailed characterization of the copy number profile of the input data. It saves BED format files of both large-structural variations and focal alterations that can be uploaded into UCSC Genome Browser.  In addition, CNAtra incorporates a visual platform that allows manual review of the CNV profile and its accessory information such as the mappability score. The detailed description of CNAtra inputs, parameters, methods, and outputs are provided in **CNAtra_User_Guide.pdf"** file.

**CNAtra output at a glance**

User can run CNAtra tool and get the CNV profile quickly using the default parameters if the BAM file is available. We provide the simulated Chr3  of CHP-212 cell line as an example. Input BAM is available as "CHP212_chr3_Artificial.bam" file under "CNAtraInput" folder. CNAtra outputs are provided in "CNAtraResults" folder (default output directory). After cloning CNAtra git repository, user can run the following commands using MATLAB to get the CNV profile.

​			*>> CNAtraDirectory = './CNAtraTool';*

​			*>> addpath(CNAtraDirectory);*

​			*>> inputFile = './CNAtraInput/CHP212_chr3_Artificial.bam';*

​			*>> CNAtraObj = CNAtra(inputFile, CNAtraDirectory);*

​			*>> CNAtraObj.RDcalculator;*

​			*>> CNAtraObj.CNVcaller;*

