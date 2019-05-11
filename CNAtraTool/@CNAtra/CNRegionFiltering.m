function [finalSegmentInfo, impRegionsPerSegment] = CNRegionFiltering(obj, targetChrIndex, signal, segmentsInfo, regionsPerSegment, refIndices, centroTeloBoundaries)
%Filtering out the significant regions around centromeres, telomeres, and black-listed areas


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ----------------------------------------- Parameters ---------------------------------------- %%
% Signal Normalization
signal = signal/obj.CNReference*2;

minImpCategory = obj.minAlterationRank; % minimum Cateogry for a region to be considered as alternation
blackGabPercentageTh = obj.maximumFalseBinsAllowed;%%Allowed region percentage to keep alteration regions
minRegionWidthForTestInBins = 15;

if(obj.removeShortRegions == 1)
	minRegionSize   = obj.resolution;
	lowerWidthLimit = minRegionWidthForTestInBins;
else
	minRegionSize   = minRegionWidthForTestInBins;% minimum region size to apply the t-test
	lowerWidthLimit = minRegionSize;
end



%%%Centromeres/Telomeres
centroStart  = centroTeloBoundaries(1,1);
centroEnd    = centroTeloBoundaries(1,2);
telo1Start   = centroTeloBoundaries(2,1);
telo1End     = centroTeloBoundaries(2,2);
telo2Start   = centroTeloBoundaries(3,1);
telo2End     = centroTeloBoundaries(3,2);



%% ------------- Regions Filtering --------------- %%
impRegionsPerSegment = containers.Map({1},{[]});
remove(impRegionsPerSegment, 1);
finalSegmentInfo = containers.Map({1},{[]});
remove(finalSegmentInfo, 1);

noSegments = length(cell2mat(keys(segmentsInfo)));
%%%%%
previousSegmentEnd = 0;
for i = 1: noSegments
	%%%%%%%%%%%%%%%%%%%%%%%%% Segments %%%%%%%%%%%%%%%%%%%%%%%%%%
	segmentInfo = segmentsInfo(i);
    refSegmentStart = refIndices(segmentInfo(1));
	refSegmentEnd   = refIndices(segmentInfo(2));
    %%
	if(refSegmentStart - previousSegmentEnd ~= 1)
		refSegmentStart = previousSegmentEnd + 1;
	end
	previousSegmentEnd = refSegmentEnd;

    %%
	refSegmentWidth = refSegmentEnd - refSegmentStart + 1;
	segmentCN    = segmentInfo(4);
	finalSegmentInfo(i) = [refSegmentStart, refSegmentEnd, refSegmentWidth, segmentCN];
	

	%%%%%%%%%%%%%%%%%%%%%%%%% Regions %%%%%%%%%%%%%%%%%%%%%%%%%%%
	regionArray = regionsPerSegment(i);
    regionArray = optimizeChangePoints(signal, segmentCN, regionArray);
    %%
    [noRegions, noParameters] = size(regionArray);
    if(noRegions > 0)
            %%---------  regions --------%%
            regionStart = regionArray(:,1);
            regionEnd   = regionArray(:,2);
            regionWidth = regionArray(:,3);
            regionRD    = regionArray(:,4);
            regionCategory = regionArray(:,5);

            %%------- Conditions --------%%
            %Condition#1: Cateogry > minimum cateogry
            impRegionsIndex = find(abs(regionCategory) >= minImpCategory);

            %Condition#2: Region includes centromeres, telomeres, gabs or black-listed area.
            thresholdBG = (1/(1-blackGabPercentageTh));
            refRegionStart = refIndices(regionStart);
            refRegionEnd   = refIndices(regionEnd);
            refRegionWidth = refRegionEnd - refRegionStart + 1;
            case1 = (refRegionWidth <= (thresholdBG*regionWidth));
            case2 = (refRegionWidth == regionWidth);
            nonCTBRegionsIndex = find((case1|case2)==1);


            %%--- Filtering regions around Black-Listed, Centromeres, and Telomeres ---------%%
            selectedRegionsIndex = intersect(impRegionsIndex, nonCTBRegionsIndex);
            impRegionsPerSegmentArray = [refRegionStart(selectedRegionsIndex), refRegionEnd(selectedRegionsIndex), refRegionWidth(selectedRegionsIndex), regionRD(selectedRegionsIndex), regionCategory(selectedRegionsIndex)];
            [noImpRegions, ~] = size(impRegionsPerSegmentArray);



            %%--- Filtering regions near Centromeres, near Telomeres, and short regions ---------%%
            filteredRegionsPerSegmentArray = [];
            for j = 1: noImpRegions
                iRegionStart = impRegionsPerSegmentArray(j,1);
                iRegionEnd   = impRegionsPerSegmentArray(j,2);
                iRegionWidth = impRegionsPerSegmentArray(j,3);
                iRegionRD    = impRegionsPerSegmentArray(j,4);
                iRegionCategory = impRegionsPerSegmentArray(j,5);

                %++ Filtering Condition ++%
                %Centromere condition
                cond1 = ((iRegionStart >= centroStart) && (iRegionEnd <= centroEnd))||((iRegionStart <= centroStart) && (iRegionEnd >= centroEnd));
                %Telomere#1 condition
                cond2 = (iRegionStart >= telo1Start)  && (iRegionEnd <= telo1End);
                %Telomere#2 condition
                cond3 = (iRegionStart >= telo2Start)  && (iRegionEnd <= telo2End);

                %Region Width
                cond4 = ((iRegionWidth <= lowerWidthLimit)) || ((iRegionWidth <= minRegionSize) && (abs(iRegionRD - segmentCN) <= 1));
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                filteringCond = cond1 || cond2 || cond3 || cond4;
                if (filteringCond == 0)
                        filteredRegionsPerSegmentArray = [filteredRegionsPerSegmentArray; impRegionsPerSegmentArray(j,:)];
                end
            end
            %%%
            impRegionsPerSegment(i) = filteredRegionsPerSegmentArray;
    else
            impRegionsPerSegment(i) = [];
    end
    %%%

end
%%%



end
%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nRegionArray = optimizeChangePoints(signal, segmentCN, regionArray)
%%%%%%%%

[noRegions, noParameters] = size(regionArray);
%
nRegionsStart = [];
nRegionsEnd   = [];
nRegionsCN    = [];

for i = 1:noRegions
    regionStart = regionArray(i,1);
    regionEnd   = regionArray(i,2);
    regionCN    = regionArray(i,4);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%% the start-point %%%%%%%%%%%%%%%%%%%%%%%%%%
    %1: left-direction
    if(regionStart > 1)
        for j = regionStart-1:-1:1
            nPointCN = signal(j);
            segmentDiff = abs(segmentCN - nPointCN);
            regionDiff  = abs(regionCN  - nPointCN);
            extendCond  = (regionDiff < segmentDiff) && (regionDiff < 1);
            if(~extendCond)
                break;
            end
        end
        nPoint = j+1;
    else
        nPoint = regionStart;
    end
    %
    %2: right-direction
    if(nPoint == regionStart)
        for j = regionStart:1:length(signal)
            nPointCN = signal(j);
            segmentDiff = abs(segmentCN - nPointCN);
            regionDiff  = abs(regionCN  - nPointCN);
            shrinkCond  = (regionDiff > segmentDiff);
            if(~shrinkCond)
                break;
            end
        end  
    end
    %
    nRegionStart = j;

    %%%%%%%%%%%%%%%%%%%%%%%%%% the end-point %%%%%%%%%%%%%%%%%%%%%%%%%%
    %1: right-direction
    if(regionEnd < length(signal))
        for j = regionEnd+1:1:length(signal)
            nPointCN = signal(j);
            segmentDiff = abs(segmentCN - nPointCN);
            regionDiff  = abs(regionCN  - nPointCN);
            extendCond  = (regionDiff < segmentDiff) && (regionDiff < 1);
            if(~extendCond)
                break;
            end
        end
        nPoint = j-1;
    else
        nPoint = regionEnd;
    end
    %
    %2: left-direction
    if(nPoint == regionEnd)
        for j = regionEnd:-1:1
            nPointCN = signal(j);
            segmentDiff = abs(segmentCN - nPointCN);
            regionDiff  = abs(regionCN  - nPointCN);
            shrinkCond  = (regionDiff > segmentDiff);
            if(~shrinkCond)
                break;
            end
        end
    end
    %
    nRegionEnd = j;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nRegionCN = median(signal(nRegionStart:nRegionEnd));

    nRegionsStart = [nRegionsStart; nRegionStart];
    nRegionsEnd   = [nRegionsEnd; nRegionEnd];
    nRegionsCN    = [nRegionsCN; nRegionCN];

end
%%%

if(~isempty(nRegionsStart))
    nRegionsWidth = nRegionsEnd - nRegionsStart+1;
    nRegionsCateogry = regionArray(:,5);
    nRegionArray = [nRegionsStart, nRegionsEnd, nRegionsWidth, nRegionsCN, nRegionsCateogry];
else
    nRegionArray = [];
end
%%%


end
%%%




