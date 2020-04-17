# CNAtra: An analytical and visualization tool for hierarchical CNA discovery of large-scale and focal copy number alterations in low-coverage cancer genomes

**CNAtra** is a MATLAB-based tool that accepts BAM/SAM files as input. It can handle both single-end and paired-end WGS reads. **CNAtra** pipeline has two modules: 1. Read depth calculator that computes, filter and normalized the RD signal,  and 2. CNV caller that identifies the copy number alterations at different length scales. **CNAtra** generates many output files providing the detailed characterization of the copy number profile of the input data. It saves BED format files of both large-scale copy number variations (LCVs) and focal alterations (FAs) that can be uploaded into UCSC Genome Browser. In addition, **CNAtra** incorporates a visual platform that allows manual review of the CNV profile and its accessory information such as the mappability score. 

**For a full description of the method and applications, please visit [CNAtra Manuscript](https://rdcu.be/b3Cki).**

**The detailed description of CNAtra software (inputs, parameters, methods, and outputs) is provided in [CNAtra user manual](CNAtra_User_Guide.pdf).**

  
## Contents
- [Download](#Download)
- [Annotations](#Annotations)
- [Parameters](#Parameters)
- [CNAtra output at a glance](#CNAtra_Example)
- [CNAtraLite](#CNAtraLite)
     
### <a name="Download"></a>Download
```bash
cd ~
git clone https://github.com/AISKhalil/CNAtra.git
cd CNAtra
```
   
### <a name="annotations"></a>Annotations  
- **CNAtra** uses the unique mappability tracks for computing the mappability scores that are used for correcting the RD signal ([Unique mappability tracks for several species](https://sites.google.com/site/anshulkundaje/projects/mappability)). This tracks includes per-base unique mappability tracks for a large range of read lengths for several key species ([Unique mappability](https://academic.oup.com/nar/article/46/20/e120/5086676)). 
- Additionally, **CNAtra** computes the GC score, for normalizing the RD signal, from the Chris Miller's pre-calculated tracks ([GC tracks](https://xfer.genome.wustl.edu/gxfer1/project/cancer-genomics/readDepth/index.html)). This includes the GC tracks for a large range of read lengths for human genome that are used for [ReadDepth software](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0016327). 
- Also, **CNAtra** uses the information of centromeres, telomeres, gap regions, and black listed regions for filtering the RD signal and CNA calls.

***N.B.***

By default, **CNAtra** uses GRCh37/hg19 human genome as the reference genome. All hg19 annotation files are included in CNAtraTool folder (`CNAtraTool/annotations`).

For other genome, user can download and build the annotation files. 

- We provided scripts for parsing the mappability and GC tracks [manual annotation](CNAtraTool/annotations/example/example.m). 
- For other annotation files (centromere.bed, telomere.bed, gap.bed, and blackListed.bed), user can make them using the tab-delimited format (chr name, start-bp, end-bp). Centromoeres' information should be provided for running **CNAtra**, while other files can be empty.

### <a name="Parameters"></a>Parameters
The main parameters of **CNAtra**.
please check [CNAtra user manual](CNAtra_User_Guide.pdf) for detailed explaination of CNAtra parameters.

    bamFile                  - the input BAM file (SAM file can be used).
   
    ploidyLevel              - the whole-genome ploidy level of the input cell ['free' (default), 'diploid', 
                             'triploid', and 'tetraploid'].
 
    minimumIBsize            - minimum segment size to be considered as Iso-copy numeric block (in kbs), 
                             default size = 1000 kb.

    resolution              - the minimum estimated-size of focal alteration based on the data coverage. 
                            User can manually tune it before before calling the ‘CNVcaller’ module.

    amplificationThreshold  - the threshold for identifying class2 focal amplfications based on the data
                            coverage. User can manually tune it before before calling the ‘CNVcaller’ module.

    deletionThreshold       - the threshold for identifying class2 focal deletions based on the data
                            coverage. User can manually tune it before before calling the ‘CNVcaller’ module.                            
### <a name="CNAtra_Example"></a>CNAtra output at a glance
User can run **CNAtra** tool and get the CNA profile quickly using the default parameters and hg19 annotation files if the BAM file is available. We provided the simulated data where we artificially incorporated LCVs and FAs in Chr3 of CHP-212 cell line data as an example for testing the tool installation. Input BAM is available as "CHP212_chr3_Artificial.bam" file under `CNAtraInput/` folder. **CNAtra** outputs are provided in `CNAtraResults/` folder (default output directory). 

Start Matlab, then edit and run the following set of commands based on your data.
```bash
% Add 'CNAtraTool' directory to Matlab search path
CNAtraDirectory = './CNAtraTool';
addpath(CNAtraDirectory);

% Define the input file
inputFile = './CNAtraInput/CHP212_chr3_Artificial.bam';

% Create CNAtra object 'CNAtraObj' with the defined parameters
CNAtraObj = CNAtra(inputFile, CNAtraDirectory);

% run 'RD calculator' module (Pipeline stage 1)
CNAtraObj.RDcalculator;

% adjust the CNAtra CNA-associated parameters before running the CNV caller module (optional)
CNAtraObj.amplificationThreshold = 0.9;
CNAtraObj.deletionThreshold = 0.9;
CNAtraObj.ploidyLevel = 'diploid';

% run 'CNV caller' module (Pipeline stage 2)
CNAtraObj.CNVcaller;

% save the CNAtraObject, so you can load it directly for further analysis.
save('CNAtraResults/CHP212_WGS.mat');

% plot the CNA tracks (e.g chr3)
chrNumber = 3;
CNAtraObj.CNAPlot('plot',chrNumber);
```
                  
### <a name="CNAtraLite"></a>CNAtraLite
**CNAtraLite** is a modified version of **CNAtra**. In **CNAtraLite**, we don't include the modules for filtering and correcting the RD signal as well as filtering of the CNA calls. Therefore, it doesn't require any annotation file (only input BAM file). This can be used for identifying the CNAs if the annotation files are not available.

Start Matlab, then edit and run the following set of commands based on your data.
```bash
% Add 'CNAtraTool' directory to Matlab search path
CNAtraDirectory = './CNAtraTool';
addpath(CNAtraDirectory);

% Define the input file
inputFile = './CNAtraInput/CHP212_chr3_Artificial.bam';

% Create CNAtraLite object 'CNAtraObj' with the defined parameters
CNAtraObj = CNAtraLite(inputFile);
  
% run 'RD calculator' module (Pipeline stage 1)
CNAtraObj.RDcalculator;

% adjust the CNAtra CNA-associated parameters before running the CNV caller module (optional)
CNAtraObj.amplificationThreshold = 0.9;
CNAtraObj.deletionThreshold = 0.9;
CNAtraObj.ploidyLevel = 'diploid';

% run 'CNV caller' module (Pipeline stage 2)
CNAtraObj.CNVcaller;

% save the CNAtraObject, so you can load it directly for further analysis.
save('CNAtraResults/CHP212_WGS.mat');

% plot the CNA tracks (e.g chr3)
chrNumber = 3;
CNAtraObj.CNAPlot('plot',chrNumber);
```
