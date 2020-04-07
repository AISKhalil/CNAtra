function [mapTracks] = mapTracks(filesDirectory, readLength, outputMat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Download unique mappability tracks of the targeted genome from https://personal.broadinstitute.org/anshul/projects/umap/
% Then, extract the .tgz file.
% Then, apply this function mapTrackForReadLength(filesDirectory, readLength) using the read length
%
binSize = 1000;
%
mapTracks = containers.Map({1},{[]});
remove(mapTracks,1);
%
for i=1:23
    targetChrIndex = i	
    %%%%%
    if(i == 23)
        j = 'X';
    else
        j = int2str(i);
    end
    
    mapFile = strcat(filesDirectory,'/chr', j,'.uint8.unique');
    tmp_uMap = fopen(mapFile,'r');
    uMapdata = fread(tmp_uMap,'*uint8');
    fclose(tmp_uMap);
    mapCond = (uMapdata > 0 & uMapdata<=readLength);
    %---------------------------------------------------------------------%
    a = mapCond;
    n = 1000;
    binnedData = arrayfun(@(k) sum(a(k:min(k+n-1,length(a)))),1:n:length(a))';
    %---------------------------------------------------------------------%
    mapTracks(i) = binnedData;
end

save(outputMat,'mapTracks');
