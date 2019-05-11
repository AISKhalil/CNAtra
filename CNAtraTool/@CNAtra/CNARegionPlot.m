function obj = CNARegionPlot(obj, saveResult, varargin)  
%Plotting the CNAtra result of a CNV region using its order in the chromosome.  

%% ------------- Plot Mode -------------- %%
nargin = length(varargin);
switch nargin,
    case 0,
        error('Error: No enough arguments.');    
    case 1,
        error('Error: No enough arguments.');    		
    case 2,
		chrNo = varargin{1};
		regionNo = varargin{2};
    otherwise
        error('Error: Too many arguments.')
end

if(strcmp(saveResult,'save')== 1)
	saveResult = 1;
else
	saveResult = 0;
end


%% ---------------- Data --------------- %%
chrData = obj.chrDictionary ;
chrNames = obj.chrNames;
segmentsInfoDic = obj.segmentsInfoDic;
regionsInfoDic  = obj.regionsInfoDic;
CNReference = obj.CNReference;
chrCentroTeloBoundaries  = obj.chrCentroTeloBoundaries;
chrBlackListedBoundaries = obj.chrBlackListedBoundaries;
chrGapBoundaries = obj.chrGapBoundaries;
chrMappabilityTracks = obj.chrMappabilityTracks;
removeLowMappabilityBins = obj.removeLowMappabilityBins;
minMappabilityThreshold = obj.minMappabilityThreshold;
binSize = obj.binSize;


%% ----------- Processing --------------%%
targetChrIndex = chrNo;
if(targetChrIndex == 23)
  targetChr = 'chrX';
else
  targetChr = strcat('chr',int2str(targetChrIndex));
end

%%%%%%%%%%%%%% Data Reading %%%%%%%%%%%%%%
orgData = chrData(targetChrIndex);
copyNumber = CNReference;
normData = orgData*2/copyNumber;
clippedDataLen = length(normData);


%%------ Segment Data  ------------%%
segmentsArray = segmentsInfoDic(targetChrIndex);
[noSegments,~] = size(segmentsArray);
segmentsStart = segmentsArray(:,2);
segmentsEnd   = segmentsArray(:,3);
segmentsCN     = segmentsArray(:,5);

segmentsCNVector = nan(clippedDataLen,1);
for j=1:noSegments
	segmentsCNVector(segmentsStart(j):segmentsEnd(j)) = segmentsCN(j);
end

%%---------- Clipped Data ----------%%
%Clip the Read-Depth data to a maximum copyNumber of 10.
maxCN = max(max(segmentsCN)+2, 10);
clippedData = min(normData, maxCN);



%%------ Regions Data  ------------%%
regionsArray = regionsInfoDic(targetChrIndex);
[noRegions,~] = size(regionsArray);
regionsStart  = regionsArray(:,2);
regionsEnd    = regionsArray(:,3);
regionsWidth   = regionsArray(:,4);
regionsCN     = regionsArray(:,5);
regionsCategory = regionsArray(:,6);

regionsCNData   = nan(clippedDataLen,1);
regionsCNVector = segmentsCNVector;
regionsCategoryVector = zeros(clippedDataLen,1);
for j=1:noRegions
	regionsCNData(regionsStart(j):regionsEnd(j)) = clippedData(regionsStart(j):regionsEnd(j));
	regionsCNVector(regionsStart(j):regionsEnd(j)) = regionsCN(j);
	regionsCategoryVector(regionsStart(j):regionsEnd(j)) = regionsCategory(j);
end


%%---- BlackListed/Gap Regions ----%%
regionsBlackGapData = nan(clippedDataLen,1);
blackListedRegions = chrBlackListedBoundaries(targetChrIndex);
gapRegions = chrGapBoundaries(targetChrIndex);

[blackRegionsNo,~] = size(blackListedRegions);
for j=1:blackRegionsNo
  regionsBlackGapData(blackListedRegions(j,1):blackListedRegions(j,2)) = clippedData(blackListedRegions(j,1):blackListedRegions(j,2));
end

[gapRegionsNo,~] = size(gapRegions);
for j=1:gapRegionsNo
  regionsBlackGapData(gapRegions(j,1):gapRegions(j,2)) = clippedData(gapRegions(j,1):gapRegions(j,2));
end

minMapValue = minMappabilityThreshold;
mapTracksCond = find(chrMappabilityTracks(targetChrIndex) <= minMapValue);
regionsBlackGapData(mapTracksCond) = clippedData(mapTracksCond);


%%------ Telomeres/Centromeres ------%%
centroTeloLocations = chrCentroTeloBoundaries(targetChrIndex);



%%%%%%%%%%%%%%%%%%%%%%%%%% Area Selecting for plot%%%%%%%%%%%%%%%%%%%%%%%%%%%
regionStart = regionsStart(regionNo);
regionEnd   = regionsEnd(regionNo);
regionWidth = regionsWidth(regionNo);

%
maxRegionCN = min(max(clippedData(regionStart:regionEnd)),maxCN);
minRegionCN = max(min(clippedData(regionStart:regionEnd)),0);
areaStart = regionStart-20*regionWidth;
areaEnd   = regionEnd + 20*regionWidth; 
plotStart = max(1, areaStart);
plotEnd   = min(areaEnd, clippedDataLen);

plotIndices = plotStart:plotEnd;

%%%%%%%%%%%%%%%%%%%%%% Ploting Copy-Number regions %%%%%%%%%%%%%%%%%%%%%%%%%%
if(saveResult == 1)
    figure(targetChrIndex);
    set(gcf, 'Position', [1 1 800 400]);
else
    figure;
    set(gcf, 'Position', get(0, 'Screensize'));
end  
hold on;
%Data
plot(plotIndices, clippedData(plotIndices),'LineStyle','none','Marker','.','Color', [0.4 0.4 0.4]);
%Amplification and deletion data
plot(plotIndices, regionsCNData(plotIndices),'b.','MarkerSize',3);
%BlackListed/gapRegions
plot(plotIndices, regionsBlackGapData(plotIndices),'k.','MarkerSize',4);
% Segments CN vector
plot(plotIndices, segmentsCNVector(plotIndices),'k','lineWidth',3);
%Amplification and deletion CN vector
plot(plotIndices, regionsCNVector(plotIndices),'r','lineWidth',2);
% Centromeres & Telomeres
centroTelo = centroTeloLocations(:);


% Region Boundary
%plot([regionStart, regionStart],[minRegionCN,maxRegionCN],'--k','LineWidth',1);
%plot([regionEnd, regionEnd],[minRegionCN,maxRegionCN],'--k','LineWidth',1); 


legend('Read-Depth','Alteration bins', 'BlackListed/Gap/Low Mapp bins', 'Segments Copy-Number', 'Alteration Copy-Number');
for j =1:length(centroTelo)
   linePt = centroTelo(j); 
   if (linePt >= plotStart && linePt <= plotEnd);
          plot([linePt,linePt],[0,5],'Color','g','LineWidth',3);   
		  legend('Read-Depth','Alteration bins', 'BlackListed/Gap/Low Mapp bins', 'Segments Copy-Number', 'Alteration Copy-Number','Centromeres-Telomeres');
   end
end    
  

ylim([-0.5,maxCN]);
ylabel('Estimated copy number');
xlabel('bin Number');
title(targetChr);
hold off;
%%%%%%%%%
dir = obj.outputDirectory;
if(exist(dir,'dir') ~= 7)
    mkdir(dir);
end
targetChr = chrNames(targetChrIndex);
%
[~,objFileName,~] = fileparts(obj.bamFile);
%
switch abs(obj.minAlterationRank)
	case 5
	    minimumClass = 'Class2';
	case 4
	    minimumClass = 'Eiffel';
	case 3 
	    minimumClass = 'Manhattan';
	case 2 
	    minimumClass = 'Class1';                                    
	otherwise
	    minimumClass = 'Neutral';    			
end
%
%
if(saveResult == 1)
        ff = strcat('-f',num2str(targetChrIndex));
        hh = strcat(obj.outputDirectory, '/', objFileName, '_', targetChr,'_minClass_', minimumClass, '_region_',int2str(regionNo),'_from_',num2str(regionStart),'Kb_to_',num2str(regionEnd),'Kb');
        print(ff,hh,'-dpng');
        savefig(hh);
        close All;
end
%
