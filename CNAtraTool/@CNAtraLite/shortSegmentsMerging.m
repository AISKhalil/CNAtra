function [finalSegmentStart, finalSegmentEnd] = shortSegmentsMerging(obj, orgSignal, orgSegmentStart, orgSegmentEnd, CNReference)
%Merging short regions with similar RD value to build the iso-copy numeric blocks. 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%------------------- Parameters ------------------%%%%%%%%%%
% maximum short segment (to define final segments).
userMinSegmentSize = obj.minimumIBsize;

% Signal Normalization
orgSignal = orgSignal/CNReference*2;



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
   if((segmentCN(i) == 0))
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




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%------ L2: Segmentation Merging (Iterative) -----%%%%%%%%%%%
SStart = L1SegmentStart;
SEnd   = L1SegmentEnd;

for i = 1:20
        SRDAdj = zeros(length(SStart),1);
        noSS = length(SStart);
        for j = 1:noSS
           segmentSample = signal(SStart(j):SEnd(j));
           SRDAdj(j) = median(segmentSample);
        end
        %%% Interval assignment for each segment
        [SCN, ~]= intervalAssignment(signal, SRDAdj, SStart, SEnd);
        [SStart, SEnd] = segmentMerging(SStart, SEnd, SRDAdj, SCN, userMinSegmentSize);
end    
%%%

[L2SegmentStart, L2SegmentEnd] = adjacentSegmentsMerging(SStart, SEnd, signal, userMinSegmentSize);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%---L3: Relatively Small Segmentation Merging (Iterative) ---%%%%
%% merging regions around CN interval boundaries
SStart = L2SegmentStart;
SEnd   = L2SegmentEnd;

for i = 1:20
        SRDAdj = zeros(length(SStart),1);
        noSS = length(SStart);
        for j = 1:noSS
           segmentSample = signal(SStart(j):SEnd(j));
           SRDAdj(j) = median(segmentSample);
        end
        %%% Interval assignment for each segment
        [SCN, ~]= intervalAssignment(signal, SRDAdj, SStart, SEnd);
        [SStart, SEnd] = segmentMerging(SStart, SEnd, SRDAdj, SCN, userMinSegmentSize);
        [SStart, SEnd] = smallSegmentMerging(signal, SStart, SEnd, userMinSegmentSize, 2);
end    

%%%
[L3SegmentStart, L3SegmentEnd] = adjacentSegmentsMerging(SStart, SEnd, signal, userMinSegmentSize);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%-- L4: Merging Adjacent Segments including the excluded regions ---%%

finalSegmentStart = L3SegmentStart;
finalSegmentEnd   = L3SegmentEnd;

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
    for i = 1:30
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






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------------------------------------------ Merging Function ----------------------------------------------------%%%
function [mergedSegmentStart, mergedSegmentEnd] = segmentMerging(segmentStart, segmentEnd, segmentsRDAdj, segmentsCN, shortSegTh)
%%%%%%%%

noSegments = length(segmentStart);
segmentLength = segmentEnd - segmentStart + 1;
shortSegmentsIndex = find(segmentLength <= shortSegTh);
longSegmentsIndex = find(segmentLength > shortSegTh);
noShortSegments = length(shortSegmentsIndex);
noLongSegments  = length(longSegmentsIndex);
totalSegmentsIndex = zeros(noSegments,1);

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
% 1) Long Segments: Merging subsequent long segments with same copy number.%

mergedSegmentsIndex = zeros(noLongSegments,1);
newOrProp = 1;%1:first segment is a new one, unconnected to previous segments.
segmentNumber = 1;
for j = 1:noLongSegments-1
   %%First Segment 
   firstSegmentIndex    = longSegmentsIndex(j); 
   firstSegmentCN       = segmentsCN(firstSegmentIndex);
   %%
   secondSegmentIndex   = longSegmentsIndex(j+1); 
   secondSegmentCN      = segmentsCN(secondSegmentIndex);
   merge = (firstSegmentCN == secondSegmentCN);
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

if(mergedSegmentsIndex(noLongSegments) == 0)
   mergedSegmentsIndex(noLongSegments) = segmentNumber;
   segmentNumber = segmentNumber+1;
end
totalSegmentsIndex(longSegmentsIndex) = mergedSegmentsIndex;


%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
% 2) Short Segments: merging short segments between long segments or create new segment %


lastLongSegmentIndex = longSegmentsIndex(end);
j=1;
while(j <= noShortSegments)
    k = shortSegmentsIndex(j);
    %%% j: segment index among short segments
    %%% k: segment index among all segments
  

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% one or group of short segments before the first long segment %%%%
    if(k == 1)
        nextlocalIndex = min(find(longSegmentsIndex> k));
        nextSegIndex = longSegmentsIndex(nextlocalIndex);
        nextSegCN = segmentsCN(nextSegIndex);
        %%%%
        firstShortSegIndex = k;
        lastShortSegIndex  = nextSegIndex-1;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %only one short segment before first long segment
        if(firstShortSegIndex==lastShortSegIndex)
          j = j+1;
          totalSegmentsIndex(firstShortSegIndex) = totalSegmentsIndex(nextSegIndex);           

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Many short segments of different CN before first long segment   
        else
      			currentSegIndices = firstShortSegIndex:lastShortSegIndex;
      			noGSeg = length(currentSegIndices);
            j = j + noGSeg;                 
      			currentGSegCN = segmentsCN(currentSegIndices);
      			currentGSegW  = segmentLength(currentSegIndices);          
      									  
      			%%
      			case1Flag = mergeOrKeepSeparateSegments1(currentGSegCN, currentGSegW, nextSegCN, shortSegTh);
      			%%
      			if(case1Flag == 0)
      				%%Make a new segment
      				totalSegmentsIndex(currentSegIndices)= segmentNumber;
      				segmentNumber = segmentNumber + 1;
      			else
      				%%Merge short segments to the next long one.
      				totalSegmentsIndex(currentSegIndices) = totalSegmentsIndex(nextSegIndex);
		      	end
            %%%
        end   
        %%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    %%%% short segment(or group of short segments) after last long segment %%%%
    elseif(k > lastLongSegmentIndex)
        prevlocalIndex = max(find(longSegmentsIndex< k));
        prevSegIndex = longSegmentsIndex(prevlocalIndex);
        prevSegCN = segmentsCN(prevSegIndex);
        %%%%
        firstShortSegIndex = k;
        lastShortSegIndex  = noSegments;
		
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %only one short segment after last long segment
        if(firstShortSegIndex==lastShortSegIndex)
		         j = j+1;
            totalSegmentsIndex(firstShortSegIndex) = totalSegmentsIndex(prevSegIndex); 
			
		    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		    %Many short segments of different CN after last long segment  
        else
            currentSegIndices = firstShortSegIndex:lastShortSegIndex;
            noGSeg = length(currentSegIndices);
            j = j + noGSeg;
			     currentGSegCN = segmentsCN(currentSegIndices);
			     currentGSegW  = segmentLength(currentSegIndices);           
      			%%			
      			case1Flag = mergeOrKeepSeparateSegments1(currentGSegCN, currentGSegW, prevSegCN, shortSegTh);			
      			%%
      			if(case1Flag == 0)
      				%%Make a new segment
      				totalSegmentsIndex(currentSegIndices)= segmentNumber;
      				segmentNumber = segmentNumber + 1;
      			else
      				%%Merge short segments to the next long one.
      				totalSegmentsIndex(currentSegIndices) = totalSegmentsIndex(prevSegIndex);
      			end
            %%%
    		end
		    %%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% short Segment (or a group of segments) between two long segments       
    else
        prevlocalIndex = max(find(longSegmentsIndex< k));
        prevSegIndex = longSegmentsIndex(prevlocalIndex);
        prevSegCN = segmentsCN(prevSegIndex);
        prevSegW = segmentLength(prevSegIndex);
        %%%%
        nextlocalIndex = min(find(longSegmentsIndex> k));
        nextSegIndex = longSegmentsIndex(nextlocalIndex);
        nextSegCN = segmentsCN(nextSegIndex);
        nextSegW = segmentLength(nextSegIndex); 
        %%%%
        firstShortSegIndex = k;
        lastShortSegIndex  = nextSegIndex-1;
		
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %only one short segment between two long segments
        if(firstShortSegIndex==lastShortSegIndex)
      			j = j+1;
      			currentSegCN = segmentsCN(firstShortSegIndex); 
      			currentSegW = segmentLength(firstShortSegIndex);

      			alterationCond = ((prevSegCN > currentSegCN) && (nextSegCN > currentSegCN)) || ((prevSegCN < currentSegCN) && (nextSegCN < currentSegCN));
      			if(totalSegmentsIndex(prevSegIndex)==totalSegmentsIndex(nextSegIndex))
        				%merge it to the surronding segments of same copy number
        				totalSegmentsIndex(firstShortSegIndex) = totalSegmentsIndex(prevSegIndex); 
      			elseif((alterationCond) && (prevSegW > nextSegW))
        				%merge to the previous segment since it is bigger than the next segment. (??????)
        				totalSegmentsIndex(firstShortSegIndex) = totalSegmentsIndex(prevSegIndex); 
      			elseif((alterationCond) && (prevSegW < nextSegW))	
        				%similar to next segment (??????)
        				totalSegmentsIndex(firstShortSegIndex) = totalSegmentsIndex(nextSegIndex); 
      			else 				
        				%Make a new segment since it is a transition between two long segments or wave-artifact (previous & next segments have different segment number) (??????)
        				totalSegmentsIndex(firstShortSegIndex)= segmentNumber;
        				segmentNumber = segmentNumber + 1;				
      			end
		  
		    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        %Many short segments between two long segments %%
        else
      			currentSegIndices = firstShortSegIndex:lastShortSegIndex;
      			noGSeg = length(currentSegIndices);       
            j = j + noGSeg;                 
            currentGSegCN = segmentsCN(currentSegIndices);
            currentGSegR  = segmentsRDAdj(currentSegIndices);
      			currentGSegW  = segmentLength(currentSegIndices);
            
			
			      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %short segments between two long segments of same segment number (same CN)
            if(totalSegmentsIndex(prevSegIndex)==totalSegmentsIndex(nextSegIndex))

        				case1Flag = mergeOrKeepSeparateSegments1(currentGSegCN, currentGSegW, prevSegCN, shortSegTh);			
        				%%
        				if(case1Flag == 0)
        					%%Make a new segment
        					totalSegmentsIndex(currentSegIndices)= segmentNumber;
        					segmentNumber = segmentNumber + 1;            
        					%%Assign a new segment number for next segment if previous and next segments have same segment number
        					totalSegmentsIndex(nextSegIndex) = segmentNumber;
        					segmentNumber = segmentNumber + 1;
        				else
        					%%Merge them to neighbors
        					totalSegmentsIndex(currentSegIndices)= totalSegmentsIndex(prevSegIndex);
        				end
      
			      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	  
            %short segments between two long segments of different segment number %%%%
            else

        				case2Flag = mergeOrKeepSeparateSegments2(currentGSegCN, currentGSegW, prevSegCN, nextSegCN, prevSegW, nextSegW, shortSegTh);	
        				%%
        				if(case2Flag == 0)
        					%%Make a new segment
        					totalSegmentsIndex(currentSegIndices)= segmentNumber;
        					segmentNumber = segmentNumber + 1;
        				elseif(case2Flag == 1)
        					%%Merge segments to the previous segment.
        					totalSegmentsIndex(currentSegIndices) = totalSegmentsIndex(prevSegIndex);
        				else
        					%%Merge segments to the next segment.
        					totalSegmentsIndex(currentSegIndices) = totalSegmentsIndex(nextSegIndex);
        				end
        				%%%
            end
	     		  %%%
        end
		    %%%
	  end
	  %%%
end  
%%%

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%3) Segment Merging: Merge based on the same segment number %%%%

lastSegment = max(totalSegmentsIndex);
mergedSegmentStart = zeros(lastSegment,1); 
mergedSegmentEnd   = zeros(lastSegment,1);
j = 1;
while(j<=lastSegment)
    segmentIndices = find(totalSegmentsIndex==j);
    segmentSL = min(segmentIndices);
    segmentEL  = max(segmentIndices);
    %%%%%%%%%
    mergedSegmentStart(j)= segmentStart(segmentSL);
    mergedSegmentEnd(j)= segmentEnd(segmentEL);    
    %%%%%%%%%
    j = j+1;
end
mergedSegmentStart = sort(mergedSegmentStart);
mergedSegmentEnd   = sort(mergedSegmentEnd);

end
%%% End-function %%%
%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%------------------ shortSegmentMergingCondition1 -----------------%%%%
function [mergeFlag]= mergeOrKeepSeparateSegments1(testSegCN, testSegW, refSegCN, shortSegTh)
%%%%%%%%
	uniqueSegCN = unique(testSegCN);
	noUniqueSegCN = length(uniqueSegCN);
	uniqueSegW = zeros(noUniqueSegCN,1);
	%%%
	for i = 1:noUniqueSegCN
		segCN = uniqueSegCN(i);
		segCNindices = find(testSegCN == segCN);
		segCNtotalWidth = sum(testSegW(segCNindices));
		uniqueSegW(i) = segCNtotalWidth;
	end
	%%%
	majorityWidth = max(segCNtotalWidth);
	majorityWidthIndex = find(segCNtotalWidth == majorityWidth);
	majoritySegCN = uniqueSegCN(majorityWidthIndex);
	%%%
	majorityCNisRefCN = (majoritySegCN == refSegCN);
	majoritySegAboveTh = (majorityWidth >= shortSegTh);
	%%%
	if(majorityCNisRefCN == 1)
		mergeFlag = 1;
	elseif(majoritySegAboveTh == 1)
		mergeFlag = 0;
	else
		mergeFlag = 1;
	end
	%%%
end
%%%
%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%------------------ shortSegmentMergingCondition2 -----------------%%%%
function [mergeFlag]= mergeOrKeepSeparateSegments2(testSegCN, testSegW, ref1SegCN, ref2SegCN, ref1SegW, ref2SegW, shortSegTh)
%%%%%%%%
	uniqueSegCN = unique(testSegCN);
	noUniqueSegCN = length(uniqueSegCN);
	uniqueSegW = zeros(noUniqueSegCN,1);
	%%%
	for i = 1:noUniqueSegCN
		segCN = uniqueSegCN(i);
		segCNindices = find(testSegCN == segCN);
		segCNtotalWidth = sum(testSegW(segCNindices));
		uniqueSegW(i) = segCNtotalWidth;
	end
	%%%
	majorityWidth = max(segCNtotalWidth);
	majorityWidthIndex = find(segCNtotalWidth == majorityWidth);
	majoritySegCN = uniqueSegCN(majorityWidthIndex);
	%%%
	alterationCond = ((ref1SegCN > majoritySegCN) && (ref2SegCN > majoritySegCN)) || ((ref1SegCN < majoritySegCN) && (ref2SegCN < majoritySegCN));
	majorityCNisRef1CN = (majoritySegCN == ref1SegCN);
	majorityCNisRef2CN = (majoritySegCN == ref2SegCN);
	majoritySegAboveTh = (majorityWidth >= shortSegTh);
	
	%%%
	if(ismember(ref1SegCN, uniqueSegCN) == 1)
		ref1SegTotalWidth = uniqueSegW(find(uniqueSegCN ==ref1SegCN));
	else
		ref1SegTotalWidth = 0;
	end
	%
	if(ismember(ref2SegCN, uniqueSegCN) == 1)
		ref2SegTotalWidth = uniqueSegW(find(uniqueSegCN ==ref2SegCN));
	else
		ref2SegTotalWidth = 0;
	end	
	%
	if(ref1SegTotalWidth > ref2SegTotalWidth)
		ref1MoreSimilar = 1;
	elseif(ref1SegTotalWidth < ref2SegTotalWidth)
		ref1MoreSimilar = 0;
	elseif(ref1SegW > ref2SegW)
		ref1MoreSimilar = 1;
	else
		ref1MoreSimilar = 0;
	end
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if(majorityCNisRef1CN == 1)
		mergeFlag = 1;
	elseif(majorityCNisRef2CN == 1)
		mergeFlag = 2;
	elseif(majoritySegAboveTh == 1)
		mergeFlag = 0;
	elseif(alterationCond == 1)
		%%Merge alteration condition (?????)
		mergeFlag = ref1MoreSimilar;
	else
		mergeFlag = 0;
	end	
	%%%
end
%%%
%%%











%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%---------------------------------------------------- Small Segment Merging -----------------------------------------%%
function [finalSegmentStart, finalSegmentEnd] = smallSegmentMerging(signal, segmentStart, segmentEnd, longSegTh, method)

segmentsWidth = segmentEnd - segmentStart;
noSegments = length(segmentStart);
segmentRDAdj = zeros(noSegments,1);

for j = 1:noSegments
   segmentSample = signal(segmentStart(j):segmentEnd(j));
   segmentRDAdj(j) = median(segmentSample);
end
[segmentsCN, ~]= intervalAssignment(signal, segmentRDAdj, segmentStart, segmentEnd);

finalSegmentStart = segmentStart(1);
finalSegmentEnd   = segmentEnd(1);
finalSegmentCN    = segmentsCN(1);
k = 1;%%Previous Index
j = 2;

while(j<noSegments)
    segmentWidth = segmentsWidth(j);
    segmentCN = segmentsCN(j);
    cSegmentRDAdj = segmentRDAdj(j);
    %%%
    prevSegmentWidth = segmentEnd(k) - segmentStart(k) + 1;
    segmentSample = signal(segmentStart(k):segmentEnd(k));
    prevSegmentRDAdj = median(segmentSample); 
    prevSegmentCN = finalSegmentCN(k);
    %%%
    nextSegmentWidth = segmentsWidth(j+1);
    nextSegmentRDAdj = segmentRDAdj(j+1);
    nextSegmentCN = segmentsCN(j+1);


    %%%%%----------------- probably an alteration region -----------------%%%%%
    cond0 = ((prevSegmentWidth > longSegTh) || (nextSegmentWidth > longSegTh));
    cond1 = ((segmentWidth/(prevSegmentWidth+nextSegmentWidth)) < 0.05); %This might be a focal alteration.
    cond2 = ((prevSegmentRDAdj > cSegmentRDAdj) && (nextSegmentRDAdj > cSegmentRDAdj))||((prevSegmentRDAdj < cSegmentRDAdj)&&(nextSegmentRDAdj < cSegmentRDAdj));%greater or smaller than its two neighbors (alteration condition)
    cond3 = (prevSegmentCN == nextSegmentCN); % similar neighbors
    %%
    condA = cond0 && cond1 && cond2 && cond3;

    %%%%%--------------- probably a false transition region --------------%%%%%
    cond4 = ((prevSegmentRDAdj > cSegmentRDAdj) && (nextSegmentRDAdj < cSegmentRDAdj))||((prevSegmentRDAdj < cSegmentRDAdj)&&(nextSegmentRDAdj > cSegmentRDAdj));%segment RD between neighbor segments' RDs
    cond5 = (segmentWidth < longSegTh) && (prevSegmentWidth > longSegTh) && (nextSegmentWidth > longSegTh);
    cond6A = ((segmentWidth/prevSegmentWidth) < 0.1) && (abs(prevSegmentRDAdj - cSegmentRDAdj) < 1);
    cond6B = ((segmentWidth/nextSegmentWidth) < 0.1) && (abs(nextSegmentRDAdj - cSegmentRDAdj) < 1);
    cond7 = abs(prevSegmentRDAdj - cSegmentRDAdj) < abs(nextSegmentRDAdj - cSegmentRDAdj);%Similar to previous segment     
    %%
    condB = cond4 && cond5 && cond6A && cond7;
    condC = cond4 && cond5 && cond6B && (~cond7);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %merge two neighbor segments
    if (condA ==1)
        finalSegmentEnd(k) = segmentEnd(j+1);
        j = j+1;

    %merge to previous segment
    elseif(condB == 1) 
        finalSegmentEnd(k) = segmentEnd(j);

    %merge to next segment
    elseif(condC == 1) 
        finalSegmentStart = [finalSegmentStart;  segmentStart(j)];
        finalSegmentEnd   = [finalSegmentEnd;  segmentEnd(j+1)];
        finalSegmentCN    = [finalSegmentCN; segmentsCN(j)];
        k = k+1;
        j = j+1;


    %keep it as a single segment   
    else 
        finalSegmentStart = [finalSegmentStart;  segmentStart(j)];
        finalSegmentEnd   = [finalSegmentEnd;  segmentEnd(j)];
        finalSegmentCN    = [finalSegmentCN; segmentsCN(j)];
        k = k+1;
    end
    %%%
    
    j = j+1;
end    

if(finalSegmentEnd(end) ~= segmentEnd(end))
        finalSegmentStart = [finalSegmentStart;  segmentStart(end)];
        finalSegmentEnd   = [finalSegmentEnd;  segmentEnd(end)];
        finalSegmentCN    = [finalSegmentCN; segmentsCN(end)];
end    

%%%
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%