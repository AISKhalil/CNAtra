function [segmentStart, segmentEnd, longSegmentStart, longSegmentEnd] = rdSavitzkySegmentation(obj, signal)
%Segmentation of the RD signal using Savitzky-Golay filter and modified-Varri techniques.

CNVsSegmentsFrequency = round(obj.resolution/0.1);%Assuming that 90% of the genome is free of focal CNVs.
minimumIBsize = obj.minimumIBsize;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%------------ Short Segments (Alteration Candidates)------------%%%%%%%%%%%%
% use Savitzky-Golay Filter and smoothing using peak envelope.
noPts1 = round(obj.resolution);%noPts1 is used as Savitzky-Golay frame length.
%
if(mod(noPts1,2)==0)
    noPts1 = noPts1+1;
end 
%
sgfSig = sgolayfilt(signal,7, noPts1);
sgfSig = sgolayfilt(sgfSig,1, noPts1);


dataLen = length(sgfSig);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%---------------------- Segmentation (1st Level) ---------------------%%
% 1) modify and apply the varri technique to find the top-ranked edges on 
% the peaks and bottoms of the signal (edge difference = peak - bottom). 
% This avoid dectection of false edges.

% a) Peaks & bottoms finding
[pks, pksLocs] = findpeaks(sgfSig);
[btsNeg, btsLocs] = findpeaks(-sgfSig);
bts = -btsNeg;

commLength = min(length(pks), length(btsNeg));
noEdges = 2*commLength-1;
if(pksLocs(1) < btsLocs(1))
    sig1 = pks(1:commLength);
    sig2 = bts(1:commLength);
    sig1Locs = pksLocs(1:commLength);
    sig2Locs = btsLocs(1:commLength);
else
    sig2 = pks(1:commLength);
    sig1 = bts(1:commLength);
    sig2Locs = pksLocs(1:commLength);
    sig1Locs = btsLocs(1:commLength);
end;
% 2) edge difference & position calculation
% to avoid outlier, we calculate edge value as the difference between subsequent peaks and bottoms. 
% edge value = diff(peak-bottom amplitude)*diff(peak- bottom Location).
% Then, we select the top edge.
index1 = 1:2:noEdges;
index2 = 2:2:noEdges;
edgeVec  = zeros(noEdges,1);
edgeLocs = zeros(noEdges,1);
vec1 = abs(sig2-sig1).*abs(sig2Locs - sig1Locs);
vec1Locs = floor((sig1Locs+sig2Locs)/2);
edgeVec(index1) = vec1;
edgeLocs(index1) = vec1Locs;
vec2 = abs(sig1(2:end)-sig2(1:end-1)).*abs(sig1Locs(2:end)-sig2Locs(1:end-1));
vec2Locs = floor((sig1Locs(2:end)+sig2Locs(1:end-1))/2);
edgeVec(index2) = vec2;
edgeLocs(index2) = vec2Locs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%----------------- Candidate Regions Selection ----------------%%
noPts2 = CNVsSegmentsFrequency*2;%noPts2 is used as frequency of segments 


% number edges = chromosome Length (in bins)/frequency of segments.
noTopEdges = dataLen/(noPts2);
topTh = floor((1-(noTopEdges/length(edgeVec)))*100);
thresholdE3 = prctile(edgeVec,topTh);
impPeaksIndexE = (edgeVec > thresholdE3);

impPLocsE = edgeLocs(impPeaksIndexE);
segmentStart = [1; impPLocsE];
segmentEnd = [impPLocsE-1;length(signal)];









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%------------------ Large Segments (Wavelets)-------------%%%%%%%%%%%%%%%%%%%%%%%
% use Savitzky-Golay Filter and smoothing using peak envelope.
noPts1 = round(minimumIBsize/2); %noPts1 is used as Savitzky-Golay frame length.
if(mod(noPts1,2)==0)
    noPts1 = noPts1+1;
end 

sgfSig = sgolayfilt(signal,7,noPts1);
sgfSig = sgolayfilt(sgfSig,1,noPts1);
dataLen = length(signal);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%--------------- Segmentation (1st Level) ------------------%%
% 1) modify and apply the varri technique to find the top-ranked edges on 
% the peaks and bottoms of the signal (edge difference = peak - bottom). 
% This avoid dectection of false edges.

% a) Peaks & bottoms finding
[pks, pksLocs] = findpeaks(sgfSig);
[btsNeg, btsLocs] = findpeaks(-sgfSig);
bts = -btsNeg;

commLength = min(length(pks), length(btsNeg));
noEdges = 2*commLength-1;
if(pksLocs(1) < btsLocs(1))
    sig1 = pks(1:commLength);
    sig2 = bts(1:commLength);
    sig1Locs = pksLocs(1:commLength);
    sig2Locs = btsLocs(1:commLength);
else
    sig2 = pks(1:commLength);
    sig1 = bts(1:commLength);
    sig2Locs = pksLocs(1:commLength);
    sig1Locs = btsLocs(1:commLength);
end;

% 2) edge difference & position calculation
% to avoid outlier, we calculate edge value as the difference between subsequent peaks and bottoms. 
% edge value = diff(peak-bottom amplitude)*diff(peak- bottom Location).
% Then, we select the top edge.
index1 = 1:2:noEdges;
index2 = 2:2:noEdges;
edgeVec  = zeros(noEdges,1);
edgeLocs = zeros(noEdges,1);
vec1 = abs(sig2-sig1).*abs(sig2Locs - sig1Locs);
vec1Locs = floor((sig1Locs+sig2Locs)/2);
edgeVec(index1) = vec1;
edgeLocs(index1) = vec1Locs;
vec2 = abs(sig1(2:end)-sig2(1:end-1)).*abs(sig1Locs(2:end)-sig2Locs(1:end-1));
vec2Locs = floor((sig1Locs(2:end)+sig2Locs(1:end-1))/2);
edgeVec(index2) = vec2;
edgeLocs(index2) = vec2Locs;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%----------------- Wavelet Candidates Selection----------------%%
noPts2 = minimumIBsize;%noPts2 is used as frequency of segments (one segment per noPts3)
noTopEdgesLS = dataLen/(noPts2);
topThLS = floor((1-(noTopEdgesLS/length(edgeVec)))*100);
thresholdE3LS = prctile(edgeVec,topThLS);
impPeaksIndexE3LS = (edgeVec > thresholdE3LS);

impPLocsELS = edgeLocs(impPeaksIndexE3LS);
longSegmentStart = [1; impPLocsELS];
longSegmentEnd = [impPLocsELS-1;length(signal)];

