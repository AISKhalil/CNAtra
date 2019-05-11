function [segmentsInfo, regionsPerSegment] = CNRegionCalling(obj, orgSignal, orgSegmentStart, orgSegmentEnd, finalSegmentStart, finalSegmentEnd, CNReference)
%identifying the significant alteration (amplified and deletion regions) using statistically t-test and the coverage-based thresholds.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -------------------- Parameters ---------------------%%
% thresholds
coverageThAmp = obj.amplificationThreshold;
coverageThDel = obj.deletionThreshold;
resolution    = obj.resolution;
minLSforStat = obj.minimumIBsize;

% iterative t-test parameters
sampleSize   = obj.sampleSize;
noIterations = obj.noIterations;
alpha = obj.alpha;

% dictionaries
segmentsInfo = containers.Map({1},{[]});
regionsPerSegment = containers.Map({1},{[]});
remove(segmentsInfo, 1);
remove(regionsPerSegment, 1);
regionsPerSegmentInfo0 = containers.Map({1},{[]});
regionsPerSegmentInfo1 = containers.Map({1},{[]});
regionsPerSegmentInfo2 = containers.Map({1},{[]});
remove(regionsPerSegmentInfo0, 1);
remove(regionsPerSegmentInfo1, 1);
remove(regionsPerSegmentInfo2, 1);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%--------------- Region Classification ----------------%%
orgSignal = orgSignal/CNReference*2;
noFinalSegments = length(finalSegmentStart);
for j =1:noFinalSegments
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Large Segment %%%%%%%%%%%%%%%%%%%%%%%%%%%
    LSegmentStart = finalSegmentStart(j);
    LSegmentEnd   = finalSegmentEnd(j);
    LSegmentWidth = LSegmentEnd - LSegmentStart +1;
    
    %%%----------Segment Parameters -----------%%%
    sample = orgSignal(LSegmentStart:LSegmentEnd);
    LSegmentRD = median(sample);
    segmentsInfo(j) = [LSegmentStart, LSegmentEnd, LSegmentWidth, LSegmentRD];


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% regions/Segment %%%%%%%%%%%%%%%%%%%%%%%%
    sSegStartNo = min(find(orgSegmentStart >= LSegmentStart));
    sSegEndNo   = max(find(orgSegmentEnd <= LSegmentEnd));
    %%
    if(sSegEndNo >= sSegStartNo)
        IBsegmentsStart = orgSegmentStart(sSegStartNo: sSegEndNo);
        IBsegmentsEnd   = orgSegmentEnd(sSegStartNo: sSegEndNo);
        if((IBsegmentsStart(1) - LSegmentStart) > resolution)
            IBsegmentsEnd   = [IBsegmentsStart(1)-1; IBsegmentsEnd];
            IBsegmentsStart = [LSegmentStart; IBsegmentsStart];
        end
        %
        if((LSegmentEnd - IBsegmentsEnd(end)) > resolution)
            IBsegmentsStart = [IBsegmentsStart; IBsegmentsEnd(end)+1];
            IBsegmentsEnd   = [IBsegmentsEnd; LSegmentEnd];
        end
        %
        noIBsegments = length(IBsegmentsStart);

        regionArray = [];%[Start, End, Width, RD, Interval-Category, P-value].
        for k = 1:noIBsegments
            regionStart = IBsegmentsStart(k);
            regionEnd   = IBsegmentsEnd(k);
            regionWidth = regionEnd - regionStart +1;
            %
            sample = orgSignal(regionStart:regionEnd);
    	    regionRD = median(sample);
        	if(LSegmentRD < 0.5 && regionRD == 0)
        		regionCategory = 1;
        		regionPValue   = 1;
        	else
        	    [regionCategory, regionPValue] = regionClassify(sample, LSegmentRD, LSegmentWidth, minLSforStat, resolution, coverageThAmp, coverageThDel, sampleSize, noIterations, alpha);
        	end
            %
            regionParameters = [regionStart, regionEnd, regionWidth, regionRD,regionCategory, regionPValue];
            regionArray = [regionArray; regionParameters];
        end
    else
        regionArray = [];%[Start, End, Width, RD, Interval-Category, P-value].
    end
    regionsPerSegmentInfo0(j) = regionArray;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%--Significant Regions Merging --%%
%% Merging regions with same class
for j = 1:noFinalSegments
    %%%Read
    segmentInfo = segmentsInfo(j);
    %%%Processing (Merging)
    LSegmentRD = segmentInfo(4);
    [mergedRegionArray] = L1RegionMerging(orgSignal, regionsPerSegmentInfo0(j), LSegmentRD, sampleSize, noIterations, alpha);
    %%%Write-Back (After Merging)
    regionsPerSegmentInfo1(j) = mergedRegionArray;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%-- Merging short neutral region between two class2 regions --%%

newMergedCategory = 5;
for j = 1:noFinalSegments
    %%%Processing (Merging)
    regionInfo = regionsPerSegmentInfo1(j);
    for i = 1:10
       [regionInfo] = L2RegionMerging(orgSignal,regionInfo, newMergedCategory);
    end
    mergedRegionArray2 = regionInfo;
    %%%Write-Back (After Merging)
    regionsPerSegmentInfo2(j) = mergedRegionArray2; 
end
%%%
regionsPerSegment = regionsPerSegmentInfo2;

%%%
%%%
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------------------------------------------- Sub-routines -----------------------------------------------------------%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------------ Category Assignment ---------------------------%%
function [regionCategory, regionPValue] = regionClassify(sample, LSegmentRD, LSegmentWidth, minLSforStat, resolution, coverageThAmp, coverageThDel, sampleSize, noIterations, alpha)
    %%%
    sampleWidth = length(sample);
    regionRD    = median(sample);

    %%% -------------------- - Coverage-threshold test -----------------------%%%
    ampRegionRD = regionRD;
    ampThreshold = coverageThAmp;
    %%%

    if(regionRD < 0.5)
        delRegionRD = regionRD;
    else
        delRegionRD = median(sample(sample~=0));
    end
    delThreshold = coverageThDel;
    % For homogenous-deletion in LSV of CN = 1       
    if(abs(LSegmentRD - delThreshold) <= 0.1)
        delThreshold = LSegmentRD - 0.1;
    end
    %%% 

    cond1 = (regionRD >= LSegmentRD);
    if(cond1 == 1)
        cond2 = abs(ampRegionRD - LSegmentRD) >= ampThreshold;
    else
        cond2 = abs(LSegmentRD  - delRegionRD) >= delThreshold;
    end

	%%% ------------------------- Statistical t-test ------------------------%%%
	if((LSegmentWidth >= minLSforStat) && (sampleWidth >= resolution))  
    	% Iteratively run t-test on a random small sample (#bins = sampleSize) from the original region to compute the P-value (the maximum one);
    	if(sampleWidth > sampleSize)
    	    samples = [];
    	    pValues = [];
    	    for j=1:noIterations;
        		iSample = randsample(sample, sampleSize);
        		[~, p] = ttest(iSample, LSegmentRD);
        		samples = [samples; iSample];
        		pValues = [pValues; p];
    	    end;
    	    %%%%%%%%%%%%%%%%%% P-Values %%%%%%%%%%%%%%%%%%%%%
    	    maxPValueIndices = find(pValues == max(pValues));
    	    maxPValueIndex = maxPValueIndices(1);
    	    maxPValue = pValues(maxPValueIndex);
    	    %%choosen P-value 
    	    regionPValue   = maxPValue;
    	else
    	    [~, regionPValue] = ttest(sample, LSegmentRD);
    	end
    	%%P-value Modification
    	regionPValue = regionPValue*(0.9*LSegmentWidth)/sampleWidth;
	else
    	regionPValue = 1;
	end
    %%%
    cond3 = regionPValue < alpha;


    %%%%%%%%%%%%%%% Category %%%%%%%%%%%%%%%
    %% Amplification (regionRD >= segmentRD)
    if(cond1 == 1)
        if((cond2 == 1))
            regionCategory = 5;
        elseif(cond3 == 1)
            regionCategory = 2;
        else
            regionCategory = 1;
        end   
    %% Deletion (regionRD < segmentRD)
    else
        %%Big difference (< threshold)
        if((cond2 == 1))
            regionCategory = -5;
        elseif(cond3 == 1)
            regionCategory = -2;
        else
            regionCategory = 1;
        end           
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------------ Level#1 Merging  ---------------------------%%
function [mergedRegionArray] = L1RegionMerging(signal, regionArray, LSegmentRD, sampleSize, noIterations, alpha)
[noRegions,noParameters] = size(regionArray);
if(noRegions > 0)
    regionStart = regionArray(:,1);
    regionEnd   = regionArray(:,2);
    regionRD    = regionArray(:,4);
    regionCategory = regionArray(:,5);


    mergedRegionIndex = zeros(noRegions,1);
    regionNumber = 1;
    mergedRegionIndex(1) = regionNumber;
    for j=2:noRegions
    		if(regionCategory(j) == regionCategory(j-1))
    			mergedRegionIndex(j) = regionNumber;
    		else	
    			regionNumber = regionNumber + 1;
    			mergedRegionIndex(j) = regionNumber;
    		end	
    end
    %%%
    lastRegion = regionNumber;
    mergedRegionArray = zeros(lastRegion,noParameters-1);%Excluding P-value

    j = 1;
    while(j<=lastRegion)
        regionIndices = find(mergedRegionIndex == j);
        regionSL = min(regionIndices);
        regionEL  = max(regionIndices);
        %%%%%%%%%
    	  LRegionStart = regionStart(regionSL);
    	  LRegionEnd   = regionEnd(regionEL);
    	  LRegionWidth = LRegionEnd-LRegionStart+1;
    	  LRegionCategory = regionCategory(regionSL);
    	  %%%%%%%%%
    	  LRegionSignal = signal(LRegionStart:LRegionEnd);
    	  LRegionRD = median(LRegionSignal);
    	  mergedRegionArray(j,:) = [LRegionStart, LRegionEnd, LRegionWidth, LRegionRD, LRegionCategory];	
        %%%%%%%%%
        j = j+1;
    end
else
    mergedRegionArray = [];
end

end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%---------------------------- Level#2 Merging  ------------------------------%%
function [mergedRegionArray] = L2RegionMerging(signal, regionArray, newCategory)
[noRegions,noParameters] = size(regionArray);
if(noRegions > 0)
    regionStart = regionArray(:,1);
    regionEnd   = regionArray(:,2);
    regionWidth = regionArray(:,3);
    regionRD    = regionArray(:,4);
    regionCategory = regionArray(:,5);

    %%Finding Core-points%%
    impCorePts = find(abs(regionCategory) == 5);
    noImpCorePts = length(impCorePts);

    if(noImpCorePts >=2)
        mergedSegmentsIndex = zeros(noRegions,1);
        mergedSegmentsIndex(1:(impCorePts(1)-1)) = [1:(impCorePts(1)-1)];%segments before first important segment
        newOrProp = 1;%1:first segment is a new one, unconnected to previous segments.
        segmentNumber = impCorePts(1);
        for j = 1:noImpCorePts-1
           %%First Segment 
           firstSegmentIndex = impCorePts(j); 
           firstSegmentSize  = regionWidth(firstSegmentIndex);
           firstSegmentSign  = sign(regionCategory(firstSegmentIndex));
           firstSegmentCN    = round(regionRD(firstSegmentIndex));

           %%Second Segment
           secondSegmentIndex = impCorePts(j+1);
           secondSegmentSize  = regionWidth(secondSegmentIndex);
           secondSegmentSign  = sign(regionCategory(secondSegmentIndex));
           secondSegmentCN     = round(regionRD(secondSegmentIndex));

           %%Segments in-between 
           inBetweenIndices = firstSegmentIndex+1: secondSegmentIndex-1;
           inBetweenCategory = regionCategory(inBetweenIndices);
           inBetweenSize = regionEnd(secondSegmentIndex-1) - regionStart(firstSegmentIndex+1) + 1;
           amplificationRejectSum = sum(inBetweenCategory < 1);
           deletionRejectSum = sum(inBetweenCategory > 1);


           %%%%%%%%%%%%%%%%%%%
           cond0 = (firstSegmentCN == secondSegmentCN);
           cond1 = (firstSegmentSign == secondSegmentSign);
           cond2 = ((inBetweenSize/(inBetweenSize+firstSegmentSize+secondSegmentSize))<= 0.1);
           %%Amplification Condition
           amplificationRejectCond = amplificationRejectSum >=1;
           cond3 = (firstSegmentSign == 1) && (amplificationRejectCond == 0);
           %%Deletion Condition
           deletionRejectCond = deletionRejectSum >=1;
           cond4 = (firstSegmentSign == -1) && (deletionRejectCond == 0);   
           %%
           mergeAmp = cond0 && cond1 && cond2 && cond3;
           mergeDel = cond0 && cond1 && cond2 && cond4;

           if(mergeAmp || mergeDel ==1)
              %%% New two segments to be merged together 
              if(newOrProp == 1)
                  mergedSegmentsIndex(firstSegmentIndex:secondSegmentIndex) = segmentNumber;
                  segmentNumber = segmentNumber + 1;
              %%% New segment to be merged with a series of previous segments    
              else
                  mergedSegmentsIndex(firstSegmentIndex+1:secondSegmentIndex) = mergedSegmentsIndex(firstSegmentIndex);
              end
              newOrProp = 0;
           else
              %%% (no merging), a new first segment            
              if(newOrProp ==1)
                  mergedSegmentsIndex(firstSegmentIndex:secondSegmentIndex-1) = segmentNumber:segmentNumber + (secondSegmentIndex-firstSegmentIndex-1);
                  segmentNumber = segmentNumber + (secondSegmentIndex-firstSegmentIndex);
              %%% (no merging), a new second segment
              else
                  mergedSegmentsIndex(firstSegmentIndex+1:secondSegmentIndex-1) = segmentNumber:segmentNumber + (secondSegmentIndex-firstSegmentIndex-2);
                  segmentNumber = segmentNumber + (secondSegmentIndex-firstSegmentIndex-1);       
              end    
              newOrProp =1;
           end   
        end

        zeroIndex = min(find(mergedSegmentsIndex == 0));
        noRegionsRemaining = noRegions - zeroIndex + 1;
        mergedSegmentsIndex(zeroIndex:end) = segmentNumber: segmentNumber+noRegionsRemaining-1;
        %%%%%%%%%%%%%%%%%%
        mergedRegionArray = [];
        lastSegment = max(mergedSegmentsIndex);
        j = 1;
        while(j<=lastSegment)
            segmentIndices = find(mergedSegmentsIndex==j);
            regionSL = min(segmentIndices);
            regionEL  = max(segmentIndices);
            %%%%%%%%%
            if(regionSL == regionEL)
                mergedRegionArray = [mergedRegionArray; regionArray(regionSL,:)];
            else
                LRegionStart = regionStart(regionSL);
                LRegionEnd   = regionEnd(regionEL);
                LRegionWidth = LRegionEnd-LRegionStart+1;
                LRegionCategory = sign(regionCategory(regionSL))*newCategory;
                %%%%%%%%%
                LRegionSignal = signal(LRegionStart:LRegionEnd);
                LRegionRD = median(LRegionSignal);
                LRegionParameters = [LRegionStart, LRegionEnd, LRegionWidth, LRegionRD, LRegionCategory];   
                mergedRegionArray = [mergedRegionArray; LRegionParameters];
            end
            j = j +1;
        end 
    else
        %No merging
        mergedRegionArray = regionArray;
    end   
else
    mergedRegionArray = [];
end 
end 
