function [regionCN] = CNAEstimator(obj, varargin)  
%Estimate copy number of a region.  

%% ------------- Plot Mode -------------- %%
nargin = length(varargin);
switch nargin,
    case 0,
        error('Error: No enough arguments.');    
    case 1,
        error('Error: No enough arguments.');    
    case 2,
        error('Error: No enough arguments.');    
	case 3,
		chrNo     = varargin{1};
        areaStart = floor(varargin{2}/obj.binSize);
		areaEnd   = ceil(varargin{3}/obj.binSize);

    otherwise
        error('Error: Too many arguments.')
end


%% ---------------- Data --------------- %%
chrData = obj.chrDictionary ;
CNReference  = obj.CNReference;
%
targetChrIndex = chrNo;
if(targetChrIndex == 23)
	targetChr = 'chrX';
else
	targetChr = strcat('chr',int2str(targetChrIndex));
end


%% ---------- CN-Estimation ----------- %%
orgData  = chrData(targetChrIndex);
normData = orgData*2/CNReference;
dataLen  = length(normData);
%
regionStart = max(1, areaStart);
regionEnd   = min(areaEnd, dataLen);
%
regionCN = median(normData(regionStart:regionEnd));


end
%%%