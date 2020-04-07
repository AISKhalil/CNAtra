function obj= binFiltering(obj)
%Filtering the RD signal to remove gap regions, low-mappability regions, black-listed regions, centromeres, and telomeres.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% chrFDictionary: RD-signal after filtering the centromeres, telomeres, and black-listed regions.
% chrFIndex: Indices of white regions (free of centromeres, telomeres, and black-listed regions).
% chrBlackCTDRegions: Indices of black regions (centromeres, telomeres, and black-listed regions).
% chrCentroTeloBoundaries: extended regions of centromeres and telomeres.
% chrBlackListedBoundaries: extended black-listed regions.





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ----------------------------- Parameters -----------------------------%%
uncertaintyDistanceToCentro = obj.uncertaintyDistanceToCentro;%number of normal bins (lower threshold < RD values < upper threshold), around the centromeres that are considered as uncertainty regions.
uncertaintyDistanceToTelo = obj.uncertaintyDistanceToTelo;%number of normal bins (lower threshold < RD values < upper threshold), around the telomeres that are considered as uncertainty regions.
keepCNVinBCT   = obj.keepCNVinBCT;   %0: remove any region around black-listed regions.
removeLowMappabilityBins = obj.removeLowMappabilityBins;
minMappabilityThreshold = obj.minMappabilityThreshold;


blackListFile = obj.blackListFile;
centromeresFile = obj.centromeresFile;
telomeresFile = obj.telomeresFile;
gapFile = obj.gapFile;
mapMatFile = obj.mapMatFile;


binSize = obj.binSize;
chrDir = obj.chrDictionary;
chrNames = obj.chrNames;
chromosomes = obj.targetChrs;


%%%%%%%%%% Filters threshold %%%%%%%%%%%%
%++Telomeres Threshold (for normal bins)
thrLT = 10;
thrUT = 90;
%++Centromeres Threshold (for normal bins)
thrLC = 10;
thrUC = 90;

minRegionSize = 3;%minimum size of gab/black region to be filtered.


%%%%%%%%%% Input Files %%%%%%%%%%%
%--------------------------------%
%%%% -- Black-list regions -- %%%%
fileID = fopen(blackListFile,'r');
blacklistRegions = textscan(fileID, '%s%f%f');
chrBlackName = blacklistRegions{1};
chrBlackStart = blacklistRegions{2};
chrBlackEnd = blacklistRegions{3};


%%%% ------- gapFile -------- %%%%
fileID = fopen(gapFile,'r');
gapRegions = textscan(fileID, '%s%f%f');
chrGapName  = gapRegions{1};
chrGapStart = gapRegions{2};
chrGapEnd   = gapRegions{3};


%%%% --- Mappability File --- %%%%
chrMappabilityTracks = obj.chrMappabilityTracks;


%%%% Centromeres & Telomeres %%%%
fileID = fopen(centromeresFile,'r');
centromeresRegions = textscan(fileID, '%s%f%f');
chrCentroName = centromeresRegions{1};
chrCentroStart = centromeresRegions{2};
chrCentroEnd = centromeresRegions{3};

%%%%
fileID = fopen(telomeresFile,'r');
telomeresRegions = textscan(fileID, '%s%f%f');
chrTeloName = telomeresRegions{1};
chrTeloStart = telomeresRegions{2};
chrTeloEnd = telomeresRegions{3};



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Filtering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------%
obj.chrFDictionary = containers.Map({1},{[]});
obj.chrFIndex = containers.Map({1},{[]});
obj.chrGapBoundaries = containers.Map({1},{[]});
obj.chrCentroTeloBoundaries = containers.Map({1},{[]});
obj.chrBlackListedBoundaries = containers.Map({1},{[]});
%%%

remove(obj.chrFDictionary,1);
remove(obj.chrFIndex,1);
remove(obj.chrGapBoundaries,1);
remove(obj.chrCentroTeloBoundaries,1);
remove(obj.chrBlackListedBoundaries,1);
%%%

for i=1:length(chromosomes)
    targetChrIndex = chromosomes(i);
    targetChr = chrNames(targetChrIndex);
    targetChrData = chrDir(targetChrIndex);
    targetChrLength = length(targetChrData);%in unit of (bins)


    %%%%%%%%%%%%%%%%%%%%%%%%%%%% Gap Bins %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---------------------------------------------------------------------%
    targetChrGapStart = [];
    targetChrGapEnd = [];
    for j=1:length(chrGapName)
        if strcmp(targetChr,chrGapName(j))== 1
            chrGapStartRegion = max(floor(chrGapStart(j)/binSize),1);
            chrGapEndRegion = ceil(chrGapEnd(j)/binSize);
            targetChrGapStart = [targetChrGapStart; chrGapStartRegion];
            targetChrGapEnd   = [targetChrGapEnd  ; chrGapEndRegion];
        end
    end
    gapBins = ones(targetChrLength,1);
    for k = 1:length(targetChrGapStart)
        if((targetChrGapEnd(k)-targetChrGapStart(k))> minRegionSize)
            gapBins(targetChrGapStart(k):targetChrGapEnd(k))= 0;
        end
    end




    %%%%%%%%%%%%%%%%%%%%%%%%% Mappability Bins%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---------------------------------------------------------------------%
    mapTrackBins = ones(targetChrLength,1);
    minMapValue = minMappabilityThreshold;
    mapTracksCond = find(chrMappabilityTracks(targetChrIndex) <= minMapValue);
    mapTrackBins(mapTracksCond) = 0; 




    %%%%%%%%%%%%%%%%%%%%%%%%% Black-Listed Bins %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---------------------------------------------------------------------%
    %Filtering the black-listed regions based on the RD of surrounding bins
    targetChrBlStart = [];
    targetChrBlEnd = [];
    for j=1:length(chrBlackName)
        if strcmp(targetChr,chrBlackName(j))== 1
            chrBlackStartRegion = max(floor(chrBlackStart(j)/binSize));
            chrBlackEndRegion = ceil(chrBlackEnd(j)/binSize);
            %%%
            targetChrBlStart = [targetChrBlStart; chrBlackStartRegion];
            targetChrBlEnd = [targetChrBlEnd; chrBlackEndRegion]; 
        end
    end    

    blackBins = ones(targetChrLength,1);
    for k = 1:length(targetChrBlStart)
        if((targetChrBlEnd(k)-targetChrBlStart(k))> minRegionSize)
            blackBins(targetChrBlStart(k):targetChrBlEnd(k))= 0;
        end
    end


    

    %%%%%%%%%%%%%% Extented Centromeres & Telomeres Bins %%%%%%%%%%%%%%%%%%
    %---------------------------------------------------------------------%
    centroTeloBins = ones(targetChrLength,1);
    %%%
    targetChrCentroStart = [];
    targetChrCentroEnd = [];
    for j=1:length(chrCentroName)
        if strcmp(targetChr,chrCentroName(j))== 1
            targetChrCentroStart = floor(chrCentroStart(j)/binSize);
            targetChrCentroEnd = ceil(chrCentroEnd(j)/binSize);
        end
    end
    %%%
    targetChrTeloStart = [];
    targetChrTeloEnd = [];
    for j=1:length(chrTeloName)
        if strcmp(targetChr,chrTeloName(j))== 1
            TeloStartTemp = floor(chrTeloStart(j)/binSize);
            if TeloStartTemp == 0
                TeloStartTemp = 1;
            end
            TeloEndTemp = ceil(chrTeloEnd(j)/binSize);
            if TeloEndTemp > targetChrLength
                TeloEndTemp = targetChrLength;%%%%
            end            
            targetChrTeloStart = [targetChrTeloStart; TeloStartTemp];
            targetChrTeloEnd = [targetChrTeloEnd; TeloEndTemp];
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (length(targetChrTeloStart)~=2)
        targetChrTeloStart = [1,targetChrLength-1];
        targetChrTeloEnd = [2,targetChrLength];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% Extend the centromeres %%%%%%%%%
    countTh = uncertaintyDistanceToCentro;
    %%%%% Centromeres Extension %%%%%
    firstBin = targetChrTeloEnd(1)+1;
    lastBin  = targetChrCentroStart;
    sample = targetChrData(firstBin:lastBin);
    nonZeroData = sample(sample ~= 0);   
    if(isempty(nonZeroData))
        targetChrCentroExStart= firstBin;
    else
        extensionTh = prctile(nonZeroData,thrLC);%%% threshold for extension
        extensionThUp = prctile(nonZeroData,thrUC);
        count = 0;
        for j = lastBin:-1:firstBin
            if(targetChrData(j)<= extensionTh ||targetChrData(j)>= extensionThUp)
                centroTeloBins(j)= 0;
            elseif (count == countTh)
                targetChrCentroExStart= j;
                break;
            else
                centroTeloBins(j)= 0;
                count = count +1;
            end
        end
        if(count ~=countTh)
           targetChrCentroExStart= firstBin;
        end
    end
    %%%%%%%%
    %%%%%%%%
    firstBin = targetChrCentroEnd;
    lastBin  = targetChrTeloStart(2)-1;
    sample = targetChrData(firstBin:lastBin);
    nonZeroData = sample(sample ~= 0);
    if(isempty(nonZeroData))
       targetChrCentroExEnd= lastBin;
    else
        extensionTh = prctile(nonZeroData,thrLC);%%% threshold for extension     
        extensionThUp = prctile(nonZeroData,thrUC);   
        count = 0;    
        for j = firstBin:1:lastBin
            if(targetChrData(j) <= extensionTh ||targetChrData(j) >= extensionThUp)
                centroTeloBins(j)= 0;
            elseif (count == countTh)
                targetChrCentroExEnd= j;
                break;
            else
                centroTeloBins(j)= 0;
                count = count +1;           
            end
        end
        if(count ~=countTh)
           targetChrCentroExEnd= lastBin;
        end
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Telomere#1 Extension %%%%%
    countTh = uncertaintyDistanceToTelo; 
    firstBin = targetChrTeloEnd(1);
    lastBin = targetChrCentroExStart-1;
    sample = targetChrData(firstBin:lastBin);
    nonZeroData = sample(sample ~= 0);
    extensionTh = prctile(nonZeroData,thrLT);%%% threshold for extension     
    extensionThUp = prctile(nonZeroData,thrUT);  
    
    targetChrTelo1ExStart = targetChrTeloStart(1);
    count = 0;       
    for j = firstBin:1:lastBin
        if(targetChrData(j) <= extensionTh ||targetChrData(j)>= extensionThUp)
            centroTeloBins(j)= 0;
        elseif (count == countTh)
            targetChrTelo1ExEnd= j;
            break;
        else
            centroTeloBins(j)= 0;
            count = count + 1;           
        end
    end
    if(count ~=countTh)
       targetChrTelo1ExEnd= lastBin;
    end
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Telomere#2 Extension %%%%%
    countTh = uncertaintyDistanceToTelo;
    firstBin = targetChrCentroExEnd+1;
    lastBin  = targetChrTeloStart(2);
    sample = targetChrData(firstBin:lastBin);
    nonZeroData = sample(sample ~= 0);
    extensionTh = prctile(nonZeroData,thrLT);%%% threshold for extension     
    extensionThUp = prctile(nonZeroData,thrUT);     
    
    count = 0;        
    for j = lastBin:-1:firstBin
        if(targetChrData(j)<= extensionTh ||targetChrData(j)>= extensionThUp)
            centroTeloBins(j)= 0;
        elseif (count == countTh)
            targetChrTelo2ExStart= j;
            break;
        else
            centroTeloBins(j)= 0;
            count = count + 1;           
        end
    end
    if(count ~=countTh)
       targetChrTelo2ExStart= firstBin;
    end
    targetChrTelo2ExEnd = targetChrTeloEnd(2); 
 
    centroTeloBins(targetChrCentroExStart:targetChrCentroExEnd)= 0;
    centroTeloBins(targetChrTelo1ExStart:targetChrTelo1ExEnd)= 0;
    centroTeloBins(targetChrTelo2ExStart:targetChrTelo2ExEnd)= 0;    




    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%% Bin Extracting %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%------------------------------------------------------------%%
    %% Blacklisted, Centromeres, and Telomeres %%
    if(keepCNVinBCT == 1)
        bCTBins = bitand(gapBins,centroTeloBins);
    else            
        bCTBins = bitand(gapBins,centroTeloBins);
        bCTBins = bitand(bCTBins,blackBins);
    end

    %% Low-Mappability %%
    if(removeLowMappabilityBins == 1)
        bCTBins = bitand(bCTBins,mapTrackBins);
    end

    whiteRegions = find(bCTBins==1);
    %blackCTRegions = find(bCTBins==0);    



   
    %%%%%%%%%%%%%%%%%%%%%%% Output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%-----------------------------------------------------------%%
    centroTeloBoundaries  = [targetChrCentroExStart, targetChrCentroExEnd; targetChrTelo1ExStart, targetChrTelo1ExEnd; targetChrTelo2ExStart, targetChrTelo2ExEnd];
    blackListedBoundaries = [targetChrBlStart, targetChrBlEnd];
    gapBoundaries         = [targetChrGapStart, targetChrGapEnd];
    
    obj.chrBlackListedBoundaries(targetChrIndex) = blackListedBoundaries;
    obj.chrGapBoundaries(targetChrIndex) = gapBoundaries;
    obj.chrCentroTeloBoundaries(targetChrIndex)  = centroTeloBoundaries;

    obj.chrFIndex(targetChrIndex) = whiteRegions;
    obj.chrFDictionary(targetChrIndex) = targetChrData(whiteRegions);
end