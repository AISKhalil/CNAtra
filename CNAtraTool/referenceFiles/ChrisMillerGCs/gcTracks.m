function [] = gcTracks(gcWindsDirectory, outputMat)
%%%%%%%%%%%%%%
% Extract gc-wind tracks at specific bin-size of different read lengths from 
% ReadDepth tool files: "https://github.com/chrisamiller/readDepth"
% & "https://xfer.genome.wustl.edu/gxfer1/project/cancer-genomics/readDepth/index.html"

%%%%%%%%%%%%%%
gcWinds = containers.Map({1},{[]});
remove(gcWinds,1);
binSize = 1000;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Reads = 100bps;
for i=1:23
    if(i ==23)
        fileName = strcat(gcWindsDirectory,'chrX.gc');
    else
        fileName = strcat(gcWindsDirectory,'chr',int2str(i),'.gc');
    end
    %
    fileName
    if exist(fileName, 'file') == 0
        gunzip(strcat(fileName,'.gz'));
    end
    %
    %
    fileID      = fopen(fileName,'r');
    GC_Data     = textscan(fileID, '%s');
    GC_Contents = GC_Data{1};
    %%%%
    a = nan(length(GC_Contents),1);
    for j = 1:length(GC_Contents)
        if(strcmp(GC_Contents(j),'NA')==0)
            a(j)= str2num(cell2mat(GC_Contents(j)));
        end
    end
 	%
    %
    n = binSize/100;
    binnedData = arrayfun(@(k) nanmean(a(k:min(k+n-1,length(a)))),1:n:length(a))';
    gcWinds(i) = binnedData;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save(outputMat,'gcWinds');

%gcTracks('Data/gcWinds100/', 'gcWinds.hg19.readLength100.Mat')
%gcTracks('Data/gcWinds76/',  'gcWinds.hg19.readLength76.Mat')
%gcTracks('Data/gcWinds50/',  'gcWinds.hg19.readLength50.Mat')
%gcTracks('Data/gcWinds36/',  'gcWinds.hg19.readLength36.Mat')
%gcTracks('Data/gcWinds27/',  'gcWinds.hg19.readLength27.Mat')
%gcTracks('Data/gcWinds200/', 'gcWinds.hg19.readLength200.Mat')
