function obj = CNAPlotInitialSegments(obj, saveResult, varargin)  
%Plotting the initial segments of a chromosome resulted before building the Iso-copy numeric blocks using our merging algorithm.  


%% ------------- Plot Mode -------------- %%
nargin = length(varargin);
switch nargin,
    case 0,
        error('Error: No enough arguments.');    
    case 1,
		% Mode:1, plot whole chromosome. 
		plotMode = 1;
		chrNo = varargin{1};
    case 2,
        % Mode:2, plot a segment.
		plotMode = 2;
		chrNo = varargin{1};
		segmentNo = varargin{2};
	case 3,
		% Mode:2, plot an area of a chromosome (in Kb).
		plotMode = 3;
		chrNo = varargin{1};
        areaStart = floor(varargin{2}/obj.binSize);
        areaEnd   = ceil(varargin{3}/obj.binSize);
    otherwise
        error('Error: Too many arguments.')
end

if(strcmp(saveResult,'save')== 1)
	saveResult = 1;
else
	saveResult = 0;
end



%% ---------------- Data --------------- %%
chrNames = obj.chrNames;
segmentsInfoDic = obj.intialSegmentsInfoDic;
regionsInfoDic  = obj.regionsInfoDic;
CNReference  = obj.CNReference;
%
chrCentroTeloBoundaries  = obj.chrCentroTeloBoundaries;
chrMappabilityTracks = obj.chrMappabilityTracks;
removeLowMappabilityBins = obj.removeLowMappabilityBins;
minMappabilityThreshold = obj.minMappabilityThreshold;
%
binSize = obj.binSize;

chromosomes = obj.targetChrs;
noChrs = length(chromosomes);



%% ----------- Processing --------------%%

targetChrIndex = chrNo;
if(targetChrIndex == 23)
	targetChr = 'chrX';
else
	targetChr = strcat('chr',int2str(targetChrIndex));
end
%%%%%%%%%%%%%%%%%%%%% Data Reading %%%%%%%%%%%%%%%%%%%	
chrLength = length(obj.chrDictionary(targetChrIndex));
%

orgData = obj.chrDictionary(targetChrIndex);
copyNumber = CNReference;
normData = orgData*2/copyNumber;
clippedDataLen = length(normData);



%%--------------- Segment Data  --------------%%
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


%%-------- Telomeres/Centromeres ---------%%
centroTeloLocations = chrCentroTeloBoundaries(targetChrIndex);



%%%%%%%%%%%%%%%%%%%%%%%%%% Area Selecting for plot%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(plotMode == 1)
    plotStart = 1;
    plotEnd   = clippedDataLen;
elseif(plotMode == 2)
    plotStart = segmentsStart(segmentNo);
    plotEnd   = segmentsEnd(segmentNo);
elseif(plotMode == 3)
    plotStart = max(1, areaStart);
    plotEnd   = min(areaEnd, clippedDataLen);
end
plotIndices = plotStart:plotEnd;



%%%%%%%%%%%%%%%%%%%%%% Ploting Copy-Number regions %%%%%%%%%%%%%%%%%%%%%%%%%%
if(saveResult == 1)
	figure(targetChrIndex);
set(gcf, 'Position', [1 1 1200 600]);
set(gcf,'Renderer','painters');
else
	figure;
set(gcf, 'Position', get(0, 'Screensize'));

end


hold on;
%Data
plot(plotIndices, clippedData(plotIndices),'LineStyle','none','Marker','.', 'MarkerSize', 4,'Color', [0.4 0.4 0.4]);
% Segments CN vector
plot(plotIndices, segmentsCNVector(plotIndices),'k','lineWidth',3);
    % Centromeres & Telomeres
    centroTelo = centroTeloLocations(:);

legend('Read-Depth', 'Initial Segments Copy-Number');   
for j =1:length(centroTelo)
   linePt = centroTelo(j); 
   if (linePt >= plotStart && linePt <= plotEnd);
          plot([linePt,linePt],[0,5],'Color','g','LineWidth',3);   
          legend('Read-Depth', 'Segments Copy-Number', 'Centromeres-Telomeres');
   end
end    
%plot(plotIndices, segmentsCNVector(plotIndices),'r','lineWidth',2);    

ylim([-0.5,maxCN]);
ylabel('Estimated copy number');
xlabel('bin Number');
title(targetChr);
hold off;
%%%%%%%%%%%%%%%%%
%%%
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
case 2 
    minimumClass = 'Class1';                                    
otherwise
    minimumClass = 'Neutral';    			
end
%
%
if(saveResult == 1)
    ff = strcat('-f',num2str(targetChrIndex));
    if(plotMode == 1)
        hh = strcat(obj.outputDirectory, '/', objFileName, '_', targetChr,'_initialSegments_minClass_', minimumClass);
    else
        hh = strcat(obj.outputDirectory, '/', objFileName, '_', targetChr,'_initialSegments_minClass_', minimumClass,'_from_',num2str(plotStart),'Kb_to_',num2str(plotEnd),'Kb');	
    end
    print(ff,hh,'-dpng');
    savefig(hh);
    close All;
end
%


end
