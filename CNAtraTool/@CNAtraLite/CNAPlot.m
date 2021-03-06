function obj = CNAPlot(obj, saveResult, varargin)  
%Plotting the CNAtra result of a chromosome, iso-copy numeric block, or genomic region.  


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
chrData = obj.chrDictionary ;
chrNames = obj.chrNames;
segmentsInfoDic = obj.segmentsInfoDic;
regionsInfoDic  = obj.regionsInfoDic;
CNReference  = obj.CNReference;
binSize = obj.binSize;


%% ----------- Processing --------------%%
targetChrIndex = chrNo;
targetChr = chrNames(targetChrIndex);

%%%%%%%%%%%%%%% Data Reading %%%%%%%%%%%%%%%
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
maxCN = max(max(segmentsCN+5),10);
clippedData = min(normData, maxCN);


regionsCNData   = nan(clippedDataLen,1);
regionsCNVector = segmentsCNVector;
regionsCategoryVector = zeros(clippedDataLen,1);
regionsCategoryVector3 = zeros(clippedDataLen,1);
regionsCategoryVector4 = zeros(clippedDataLen,1);
regionsCategoryVector5 = zeros(clippedDataLen,1);

%%------ Regions Data  ------------%%
regionsArray = regionsInfoDic(targetChrIndex);
[noRegions,~] = size(regionsArray);
if(noRegions > 0)
	regionsStart  = regionsArray(:,2);
	regionsEnd    = regionsArray(:,3);
	regionsCN     = regionsArray(:,5);
	regionsCategory = regionsArray(:,6);	

	for j=1:noRegions
		regionsCNData(regionsStart(j):regionsEnd(j)) = clippedData(regionsStart(j):regionsEnd(j));
		regionsCNVector(regionsStart(j):regionsEnd(j)) = regionsCN(j);
		regionsCategoryVector(regionsStart(j):regionsEnd(j)) = regionsCategory(j);
		%%
		if(abs(regionsCategory(j))==3)
			regionsCategoryVector3(regionsStart(j):regionsEnd(j)) = regionsCategory(j);
		end
		if(abs(regionsCategory(j))==4)
			regionsCategoryVector4(regionsStart(j):regionsEnd(j)) = regionsCategory(j);
		end
		if(abs(regionsCategory(j))==5)
			regionsCategoryVector5(regionsStart(j):regionsEnd(j)) = regionsCategory(j);
		end
	end
end

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
%
hold on;
%
plot(plotIndices, clippedData(plotIndices),'LineStyle','none','Marker','.', 'MarkerSize', 4,'Color', [0.4 0.4 0.4]);
%Amplification and deletion data
plot(plotIndices, regionsCNData(plotIndices),'b.','MarkerSize',3);
%Amplification and deletion CN vector
plot(plotIndices, regionsCNVector(plotIndices),'r','lineWidth',2);
% Segments CN vector
plot(plotIndices, segmentsCNVector(plotIndices),'k','lineWidth',3);

legend('Read-Depth','Alteration bins', 'Alteration Copy-Number', 'Segments Copy-Number');   

ylim([-0.5,maxCN]);
ylabel('Copy number');
xlabel('Bin Number');
title(targetChr);
hold off;
%%%%%%%%%%%%%%%%%
 
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
        hh = strcat(obj.outputDirectory, '/', objFileName, '_', targetChr,'_minClass_', minimumClass);
    else
        hh = strcat(obj.outputDirectory, '/', objFileName, '_', targetChr,'_minClass_', minimumClass,'_from_',num2str(plotStart),'Kb_to_',num2str(plotEnd),'Kb');	
    end
    print(ff,hh,'-dpng');
    savefig(hh);
    close All;
end


end
