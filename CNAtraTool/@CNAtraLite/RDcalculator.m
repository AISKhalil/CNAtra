function obj = RDcalculator(obj)
%Computing the read depth signal from the input reads.


%%%%%%%%%%%%%%%%
% input: SAM/BAM 
% output: dictionary of reads/bins for each chromosome, chromosome lengths, and chromosome names.

bamFile = obj.bamFile;
%%%%% Bam binning
[~, ~, ext] = fileparts(bamFile);
if strcmp(ext,'.bam') || strcmp(ext,'.sam')
    bamBinning(obj, bamFile);
else
    error('RDcalculator:FileExtension', ...
        'File should be in SAM/BAM Format');   
end
%   
obj.chrDictionary = obj.chrRawDictionary;

%%%%% Analysis parameters
pipelineParameters(obj);

