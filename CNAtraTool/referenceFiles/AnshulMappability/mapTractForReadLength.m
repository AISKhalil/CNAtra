function [mapTracks] = mapTrackForReadLength(filesDirectory, readLength)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
binSize = 1000;

mapTracks = containers.Map({1},{[]});
remove(mapTracks,1);


for i=1:23
    targetChrIndex = i	
    %%%%%
    if(i == 23)
        j = 'X';
    else
        j = int2str(i);
    end
    
    mapFile = strcat(filesDirectory,'chr', j,'.uint8.unique');
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


save(strcat('mapTracks.hg19.',int2str(readLength),'.mat'),'mapTracks');
