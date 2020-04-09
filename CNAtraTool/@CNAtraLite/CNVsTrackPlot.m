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
chrData = obj.chrDictionary ;
segmentsInfoDic = obj.segmentsInfoDic;
regionsInfoDic  = obj.regionsInfoDic;
CNReference = obj.CNReference;
binSize = obj.binSize;
%
chromosomes = obj.targetChrs;
noChrs = length(chromosomes);




%% ----------- Processing --------------%%
targetChrIndex = chrNo;
targetChr = chrNames(targetChrIndex);

%%%%%%%%%%%%%%%%%%%%% Data Reading %%%%%%%
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
		segmentsCNVector(regionsStart(j):regionsEnd(j)) = regionsCN(j);
	end
end


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
% Segments CN vector
plot(plotIndices, segmentsCNVector(plotIndices),'k','lineWidth',4);
legend('Read depth', 'Copy number');   

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