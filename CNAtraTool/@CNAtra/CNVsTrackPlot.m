function obj = CNVsTrackPlot(obj, saveResult, varargin)  
%Plotting the CNAtra CNV track (contructed by merging the iso-copy numeric blocks and focal amplifications/deletions) of a chromosome.  


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
		% Mode:3, plot an area of a chromosome (in Kb).
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
segmentsInfoDic = obj.segmentsInfoDic;
regionsInfoDic  = obj.regionsInfoDic;
CNReference = obj.CNReference;
%
chrCentroTeloBoundaries  = obj.chrCentroTeloBoundaries;
chrMappabilityTracks = obj.chrMappabilityTracks;
removeLowMappabilityBins = obj.removeLowMappabilityBins;
minMappabilityThreshold = obj.minMappabilityThreshold;
%
binSize = obj.binSize;
%
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
chrFDictionary = obj.chrFDictionary(targetChrIndex);
chrFIndex = obj.chrFIndex(targetChrIndex);
chrLength = length(obj.chrDictionary(targetChrIndex));
%
orgData = zeros(chrLength,1);
orgData(chrFIndex) = chrFDictionary;
%
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


%%----- Clipped Data -------%%
%Clip the Read-Depth data to a maximum copyNumber of 10.
maxCN = max(max(segmentsCN+2),10);
clippedData = min(normData, maxCN);

%%------- Regions Data  ----------%%	
regionsCNData   = nan(clippedDataLen,1);
regionsCNVector = segmentsCNVector;	
%
regionsArray = regionsInfoDic(targetChrIndex);
[noRegions,~] = size(regionsArray);
%
if(noRegions > 0)	
	regionsStart  = regionsArray(:,2);
	regionsEnd    = regionsArray(:,3);
	regionsCN     = regionsArray(:,5);
	regionsCategory = regionsArray(:,6);

	for j=1:noRegions
		regionsCNData(regionsStart(j):regionsEnd(j)) = clippedData(regionsStart(j):regionsEnd(j));
		regionsCNVector(regionsStart(j):regionsEnd(j)) = regionsCN(j);
		%
		if(abs(regionsCategory(j))==3)
			segmentsCNVector(regionsStart(j):regionsEnd(j)) = regionsCN(j);
		end
		%
		if(abs(regionsCategory(j))==4)
			segmentsCNVector(regionsStart(j):regionsEnd(j)) = regionsCN(j);
		end
		%
		if(abs(regionsCategory(j))==5)
			segmentsCNVector(regionsStart(j):regionsEnd(j)) = regionsCN(j);
		end
	end
end


%%------------- Filtered-Regions ---------%%
regionsBlackGapData = nan(clippedDataLen,1);
chrBinsNumbers = 1:clippedDataLen;
blackBinsNumbers = setdiff(chrBinsNumbers, chrFIndex);
regionsBlackGapData(blackBinsNumbers) = clippedData(blackBinsNumbers);


%%------ Telomeres/Centromeres ------%%
centroTeloLocations = chrCentroTeloBoundaries(targetChrIndex);



%%%--------- Area Selecting for plot ---------%%%
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



%%%------ Ploting Copy-Number regions -------%%%
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
plot(plotIndices, clippedData(plotIndices),'LineStyle','none','Marker','.','Color', [0.4 0.4 0.4]);
%Gap data
plot(plotIndices, regionsBlackGapData(plotIndices),'k.','MarkerSize',4);
% Segments CN vector
plot(plotIndices, segmentsCNVector(plotIndices),'k','lineWidth',4);
% Centromeres & Telomeres
centroTelo = centroTeloLocations(:);

legend('Read-Depth', 'Filtered-bins', 'Segments Copy-Number');   
for j =1:length(centroTelo)
linePt = centroTelo(j); 
   if (linePt >= plotStart && linePt <= plotEnd);
      plot([linePt,linePt],[0,5],'Color','g','LineWidth',3);   
      legend('Read-Depth', 'Filtered-bins', 'Segments Copy-Number', 'Centromeres-Telomeres');
   end
end    


ylim([-0.5,maxCN]);
ylabel('Estimated copy number');
xlabel('bin Number');
title(targetChr);
hold off;


%%%------ saving output-figure --------%%%
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
        hh = strcat(obj.outputDirectory, '/', objFileName, '_', targetChr,'_CNVsTrack_minClass_', minimumClass);
    else
        hh = strcat(obj.outputDirectory, '/', objFileName, '_', targetChr,'_CNVsTrack_minClass_', minimumClass,'_from_',num2str(plotStart),'Kb_to_',num2str(plotEnd),'Kb');	
    end
    print(ff,hh,'-dpng');
    savefig(hh);
    close All;
end
%


end
%%%