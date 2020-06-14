CNAtraDirectory = '../CNAtra/CNAtraTool';
addpath(CNAtraDirectory);

% Define the input file
inputFile = './CNAtraInput/CHP212_chr3_Artificial.bam';

% Create CNAtra object 'CNAtraObj' with the defined parameters
CNAtraObj = CNAtra(inputFile, CNAtraDirectory);
CNAtraObj.outputDirectory = './CNAtraResults';

% run 'RD calculator' module (Pipeline stage 1)
CNAtraObj.RDcalculator;

% adjust the CNAtra CNA-associated parameters before running the CNV caller module (optional)
CNAtraObj.amplificationThreshold = 0.9;
CNAtraObj.deletionThreshold = 0.9;
CNAtraObj.ploidyLevel = 'diploid';

% run 'CNV caller' module (Pipeline stage 2)
CNAtraObj.CNVcaller;

% save the CNAtraObject, so you can load it directly for further analysis.
save('./CNAtraResults/CHP212_WGS.mat');

% save the CNA tracks (e.g chr3)
chrNumber = 3;
CNAtraObj.CNAPlot('save',chrNumber);
CNAtraObj.CNVsTrackPlot('save',chrNumber);

% save the genome-wide tracks
binSize = 30000;
CNAtraObj.plotGenome('save', binSize)
CNAtraObj.plotGenomeHistogram('save', binSize)
CNAtraObj.ploidyTest('save')