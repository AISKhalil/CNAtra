function [segmentsInfo] = neighborSegmentsMerging(obj, targetChrIndex, orgSignal, orgSegmentStart, orgSegmentEnd, CNReference)
%Merging adjacent short regions with similar RD value to build the CNV tracks. 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%------------------- Parameters ------------------%%%%%%%%%%
% maximum short segment (to define final segments).
userMinSegmentSize = obj.minimumIBsize;

% Signal Normalization
orgSignal = orgSignal/CNReference*2;

% Merged segments Info
segmentsInfo = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%- L1: Merging Adjacent Segments (Excluding Homogenous deletion)-%% 

noSeg = length(orgSegmentStart);
RDAdj = zeros(length(orgSegmentStart),1);
for j = 1:noSeg
   segmentSample = orgSignal(orgSegmentStart(j):orgSegmentEnd(j));
   RDAdj(j) = median(segmentSample);
end
%

[segmentCN, ~] = intervalAssignment(orgSignal, RDAdj, orgSegmentStart, orgSegmentEnd);
segmentWidth   = orgSegmentEnd - orgSegmentStart +1;
nonHDelVec     = ones(length(orgSignal),1);
segmentStart   = [];
segmentEnd     = [];
segmentRDdd    = [];
%
offset = 0;
for i = 1:noSeg
   if( (segmentCN(i) == 0))
    offset = offset + (orgSegmentEnd(i)-orgSegmentStart(i)+1);
    nonHDelVec(orgSegmentStart(i):orgSegmentEnd(i)) = 0; 
   else
    segmentStart = [segmentStart; orgSegmentStart(i)-offset];
    segmentEnd   = [segmentEnd; orgSegmentEnd(i)-offset];
    segmentRDdd  = [segmentRDdd; RDAdj(i)];
   end     
end
%
nonHDelIndex = find(nonHDelVec==1);
signal = orgSignal(nonHDelIndex);

[L1SegmentStart, L1SegmentEnd] = adjacentSegmentsMerging(segmentStart, segmentEnd, signal, userMinSegmentSize);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%-- L2: Merging Adjacent Segments including the excluded regions ---%%

finalSegmentStart = L1SegmentStart;
finalSegmentEnd   = L1SegmentEnd;

finalSegmentStart = nonHDelIndex(finalSegmentStart);
finalSegmentEnd   = nonHDelIndex(finalSegmentEnd);
noFinalSegments = length(finalSegmentStart);
for j =1:noFinalSegments-1
   %j; 
   firstSegmentEnd   = finalSegmentEnd(j); 
   firstSegmentRD = median(orgSignal(finalSegmentStart(j):finalSegmentEnd(j)));
   %j+1; 
   secondSegmentStart   = finalSegmentStart(j+1); 
   secondSegmentRD = median(orgSignal(finalSegmentStart(j+1):finalSegmentEnd(j+1)));   
   %
   startEndDiff = secondSegmentStart - firstSegmentEnd;
   firstIsLower = (firstSegmentRD <= secondSegmentRD);
    
   if(startEndDiff > 1)
      if(firstIsLower)%Extend first segment to include deletion area
          finalSegmentEnd(j) = finalSegmentStart(j+1)-1;
      else %Extend second segment to include deletion area
          finalSegmentStart(j+1) = finalSegmentEnd(j)+1;
      end    
   end    
end

finalSegmentStart(1)=1;
finalSegmentEnd(end)= length(orgSignal); 

[finalSegmentStart, finalSegmentEnd] = adjacentSegmentsMerging(finalSegmentStart, finalSegmentEnd, orgSignal, userMinSegmentSize);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
refIndices = obj.chrFIndex(targetChrIndex);
%%
noFinalSegments = length(finalSegmentStart);
previousSegmentEnd = 0;
for j =1:noFinalSegments
    LSegmentStart = finalSegmentStart(j);
    LSegmentEnd   = finalSegmentEnd(j);
    LSegmentWidth = LSegmentEnd - LSegmentStart +1;
    sample = orgSignal(LSegmentStart:LSegmentEnd);
    LSegmentRD = median(sample);
    %%%
    chrSegmentInfo = [j, LSegmentStart, LSegmentEnd, LSegmentWidth, LSegmentRD];
    chrSegmentInfo(:,2) = refIndices(chrSegmentInfo(:,2));
    chrSegmentInfo(:,3) = refIndices(chrSegmentInfo(:,3));
    %%%
    if(chrSegmentInfo(:,2) - previousSegmentEnd ~= 1)
  chrSegmentInfo(:,2) = previousSegmentEnd + 1;
    end
    previousSegmentEnd = chrSegmentInfo(:,3);
    %%%
    segmentsInfo = [segmentsInfo; chrSegmentInfo];
end






end
%%%
%%%
%%%






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%_____________________________________________________ Sub-rountines ________________________________________________________%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%------------------------ Adjacent Segment Merging ------------------------------------%%%%%%%%%
function [SStart, SEnd] = adjacentSegmentsMerging(segmentStart, segmentEnd, signal, userMinSegmentSize)
%%%%%%%%
    SStart = segmentStart;
    SEnd = segmentEnd;
    %
    for i = 1:50
        [SStart, SEnd] = segmentRefinement(signal, SStart, SEnd);
    %%%%
        noMergedSegments = length(SStart);
        mergedSegmentRDAdj = zeros(noMergedSegments,1);
        for j = 1:noMergedSegments
            segmentSample = signal(SStart(j):SEnd(j));
            mergedSegmentRDAdj(j) = median(segmentSample);
        end
        [mergedSegmentCN, ~] = intervalAssignment(signal, mergedSegmentRDAdj, SStart, SEnd);
        [SStart, SEnd] = segmentEdgeRefinement(SStart, SEnd, mergedSegmentRDAdj, mergedSegmentCN, userMinSegmentSize);
    %%%%
        [SStart, SEnd] = segmentRefinement(signal, SStart, SEnd);
    end
end
%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%------------- Segment Refinement -----------------%%%%%%%%%
function [nSStart, nSEnd] = segmentRefinement(signal, SStart, SEnd)
j = 1;
while(j <length(SStart))
    %%%
    segmentSample1 = signal(SStart(j):SEnd(j));
    seg1Avg = median(segmentSample1);
    [seg1CN, ~] = intervalAssignment(signal, seg1Avg, SStart(j), SEnd(j));
    %%%
    segmentSample2 = signal(SStart(j+1):SEnd(j+1));
    seg2Avg = median(segmentSample2);
    [seg2CN, ~] = intervalAssignment(signal, seg2Avg, SStart(j+1), SEnd(j+1));
        
    if(seg1CN == seg2CN)
        SEnd(j)= SEnd(j+1);
        SStart(j+1) = [];
        SEnd(j+1) = [];
        j = j-1;
    end
    j = j+1;
end
nSStart = SStart;
nSEnd   = SEnd;
end
%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%------------------ Edge Refinement -----------------%%%%%%%%
function [mergedSegmentStart, mergedSegmentEnd] = segmentEdgeRefinement(segmentStart, segmentEnd, segmentsRDAdj, segmentsCN, shortSegTh)

segmentLength = segmentEnd - segmentStart+1;
noSegments = length(segmentStart);
mergedSegmentsIndex = zeros(noSegments,1);
newOrProp = 1;%1:first segment is a new one, unconnected to previous segments.
segmentNumber = 1;
for j = 1:noSegments-1

   %%First Segment 
   firstSegmentIndex = j; 
   firstSegmentCN  = segmentsCN(firstSegmentIndex);
   firstSegmentRD = segmentsRDAdj(firstSegmentIndex);
   %%
   secondSegmentIndex = j+1; 
   secondSegmentRD = segmentsRDAdj(secondSegmentIndex);
   secondSegmentCN  = segmentsCN(secondSegmentIndex);
   
   %%Conditions
   cond1 = (firstSegmentCN ~= secondSegmentCN);% Alternate-around-edge;
   crossMergingThreshold = 0.5;
   cond2 = abs(secondSegmentRD - firstSegmentRD) < crossMergingThreshold;
   merge = cond1 && cond2;
   
   if(merge ==1)
      %%% New two segments to be merged together 
      if(newOrProp == 1)
          mergedSegmentsIndex(j:j+1) = segmentNumber;
          segmentNumber = segmentNumber + 1;
      %%% New segment to be merged with a series of previous segments    
      else
          mergedSegmentsIndex(j+1) = mergedSegmentsIndex(j);
      end
      newOrProp = 0;
   else
      %%% New two segments,with different medians (no merging).            
      if(newOrProp ==1)
          mergedSegmentsIndex(j) = segmentNumber;
          segmentNumber = segmentNumber + 1;
      end    
      newOrProp =1;
   end   
end

%Last Segment
if(mergedSegmentsIndex(noSegments) == 0)
   mergedSegmentsIndex(noSegments) = segmentNumber;
end
%%%
lastSegment = max(mergedSegmentsIndex);
mergedSegmentStart = zeros(lastSegment,1); 
mergedSegmentEnd   = zeros(lastSegment,1);
j = 1;
while(j<=lastSegment)
    segmentIndices = find(mergedSegmentsIndex==j);
    segmentSL = min(segmentIndices);
    segmentEL = max(segmentIndices);
    %%%%%%%%%
    mergedSegmentStart(j)= segmentStart(segmentSL);
    mergedSegmentEnd(j)= segmentEnd(segmentEL);    
    %%%%%%%%%
    j = j+1;
end
mergedSegmentStart = sort(mergedSegmentStart);
mergedSegmentEnd   = sort(mergedSegmentEnd);
end
%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%---------------------- Interval Assignment ----------------------%%%%
function [segmentsCN, CNVector]= intervalAssignment(signal,segmentsRDAdj, segmentStart, segmentEnd) 
%%%%%%%%
      i1 = [0.5:1:99]';
      i2 = [0;i1];
      i3 = [i1;99.5];
      CNIntervels = [i2,i3]; 
      %%
      [noCNIntervals,~]=size(CNIntervels);
      noSegments = length(segmentsRDAdj);
      CNVector = zeros(length(signal),1);
      segmentsCN = zeros(noSegments,1);
      for j = 1:noSegments
        segmentRD = segmentsRDAdj(j);
        for k = 1:noCNIntervals
          CNStart = CNIntervels(k,1);
          CNEnd   = CNIntervels(k,2);
          intervalIn = (segmentRD >CNStart && segmentRD<= CNEnd);
          if(intervalIn ==1)
             segmentsCN(j) = k-1;
             break;
          end    
        end 
        CNVector(segmentStart(j):segmentEnd(j))= segmentsCN(j); 
      end
end
%%%




