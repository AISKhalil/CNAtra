% This script gives an example to download and create GC and mappability tracks
% for hg18 and read length of 36.


% 1) From "https://xfer.genome.wustl.edu/gxfer1/project/cancer-genomics/readDepth/index.html",
% - download the gc-wind file (100-bp bins) for a read length of 36 bp.
% - extract it to "./input" directory.
gcWindsDirectory = './input/gcWinds.readLength36.hg18/gcWinds';


% 2) From "https://personal.broadinstitute.org/anshul/projects/umap/", 
% - download the unique mappability tracks of the targeted genome and the read length.
% - extract it (./input/globalmap_k32tok42.tgz) to the "./input" directory.
mappDirectory = './input/globalmap_k32tok42';

% 3) Run the gcTracks.m and map.tracks to parsing these files 
readLength = 36;
gcTracksMat  = strcat('./annotatedFiles/gcWinds.',int2str(readLength),'.mat');
mapTracksMat = strcat('./annotatedFiles/mapTracks.',int2str(readLength),'.mat');
%
gcTracks(gcWindsDirectory, gcTracksMat)
mapTracks(mappDirectory, readLength, mapTracksMat)