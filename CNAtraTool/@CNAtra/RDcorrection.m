function obj = RDcorrection(obj)
% GC-contents and Mappability correction.


%%%%-----------------  Input Data ---------------%%%%
id = 'bioinfo:saminfo:InvalidTagField';
warning('off',id);
inMemoryOpt = 'false';
MAPQ_Score = obj.MAPQ_Score;
binSize    = obj.binSize;
chrNames   = obj.chrNames;
chrLengths = obj.chrLengths;

obj.chrDictionary = containers.Map({1},{[]});
remove(obj.chrDictionary,1);

obj.chrMappabilityTracks = containers.Map({1},{[]});
remove(obj.chrMappabilityTracks,1);


chrLengths = [];
for i = 1:length(obj.chrLengths.keys())
    chrLengths = [chrLengths, obj.chrLengths(obj.targetChrs(i))];
end
binsCounts = ceil(chrLengths/binSize);
%
binsStop = arrayfun(@(k) sum(binsCounts(1:k)),1:1:length(binsCounts));
binsStart = [1,binsStop(1:length(binsCounts)-1)+1];

readsPerBins = zeros(binsStop(end),1);
mappabilityScores = zeros(binsStop(end),1);
gcContents   = nan(binsStop(end),1);


%%%% --- Mappability File --- %%%%
mapMatFile = obj.mapMatFile;
h1 = load(mapMatFile);
for i = 1:1:length(chrNames)
    %%%-Chromosome Subset-%%
    targetChrIndex = obj.targetChrs(i);
    targetChr = chrNames(targetChrIndex);
    obj.chrMappabilityTracks(targetChrIndex) = h1.mapTracks(targetChrIndex)/binSize;
    startLocation = binsStart(i);
    stopLocation  = binsStop(i);    
    mappabilityScores(startLocation:stopLocation) = obj.chrMappabilityTracks(targetChrIndex);
end
clear h1;
if(obj.minMappabilityThreshold == 1)
    minMappabilityThreshold = prctile(mappabilityScores,10);
    obj.minMappabilityThreshold = minMappabilityThreshold;
end
minMappabilityThreshold = obj.minMappabilityThreshold;






%%%% -- Loading GC-contents -- %%%%
gcWindsMatFile = obj.gcWindsMatFile;
h1 = load(gcWindsMatFile);
obj.chrGCTracks = h1.gcWinds;
clear h1;
gcCorrectionMethod = obj.gcCorrectionMethod;


%%%%------------------- Correction -------------%%%%


%%%%--- Mappability Correction ---%%%%
for i = 1:1:length(chrNames)

    %%%-Chromosome Subset-%%
    targetChrIndex = obj.targetChrs(i);
    targetChr = chrNames(targetChrIndex);
    %
    binnedData = obj.chrRawDictionary(targetChrIndex);
    mappabilityTracks = obj.chrMappabilityTracks(targetChrIndex);
    MappabilityIndices = (mappabilityTracks > minMappabilityThreshold);

    %% Mappability Normalization %%
    binnedData(MappabilityIndices) = binnedData(MappabilityIndices) ./ mappabilityTracks(MappabilityIndices);
    obj.chrDictionary(targetChrIndex) = binnedData;

    % Genome-wise Data
    startLocation = binsStart(i);
    stopLocation  = binsStop(i);
    readsPerBins(startLocation:stopLocation) = binnedData;
    clear MappabilityIndices mappabilityTracks binnedData;
    
    % Genome-wise GC-contents
    if(gcCorrectionMethod == 1)
        gcContents(startLocation:stopLocation) = obj.chrGCTracks(targetChrIndex);
        if(targetChrIndex == 1)
            disp('GC-correction using ChrisMiller tracks');
        end    
    elseif(gcCorrectionMethod == 2)
        gcContents(startLocation:stopLocation) = obj.chrInputGCContents(targetChrIndex);
        if(targetChrIndex == 1)
            disp('GC-correction using data GC-contents');
        end          
    else
        if(targetChrIndex == 1)
            disp('No GC-correction');
        end        
    end
end



%%%%--- GC-Correction ---%%%%
genomeMappabilityIndices = (mappabilityScores > minMappabilityThreshold);

if(gcCorrectionMethod == 1 || gcCorrectionMethod == 2)
    reads = readsPerBins;
    readsMod = reads;
    gc = gcContents;
    valueIndex = (~ isnan(gc)) & genomeMappabilityIndices;
    rdGlobal = mean(reads(valueIndex));
    %%%
    rangeBoundaries = 0.005:0.005:0.995;
    noIntervals = length(rangeBoundaries)-1;
    for j = 1:noIntervals
       intervalIndex   = ( rangeBoundaries(j) <= gc & gc < rangeBoundaries(j+1)) & (genomeMappabilityIndices);
       if(sum(intervalIndex)>1)
           intervalData    = reads(intervalIndex);
           intervalRDMean  = mean(intervalData);
           correctionRatio = rdGlobal/intervalRDMean;
           if(intervalRDMean == 0)
                intervalDataMod = intervalData;
           else
                intervalDataMod = intervalData*correctionRatio;
           end
           readsMod(intervalIndex) = intervalDataMod;
       end 
    end 
    %%%
    for i = 1:1:length(chrNames)
        obj.chrDictionary(obj.targetChrs(i)) = readsMod(binsStart(i):binsStop(i)); 
    end    
    clear reads readsMod
else
    %No-Correction
    for i = 1:1:length(chrNames)
        obj.chrDictionary(obj.targetChrs(i)) = readsPerBins(binsStart(i):binsStop(i)); 
    end     
end
%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
