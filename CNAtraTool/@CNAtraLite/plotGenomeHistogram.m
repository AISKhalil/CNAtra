function obj = plotGenomeHistogram(obj, saveResult, nBinSize)  
%Plotting the genome RD frequency distribution.


%% -------------- Figure setup ----------- %%
if(strcmp(saveResult,'save')== 1)
    saveResult = 1;
else
    saveResult = 0;
end
%
if(saveResult == 1)
    figure(102);
    set(gcf, 'Position', [1 1 1200 600]);
    set(gcf,'Renderer','painters');
else
    figure;
    set(gcf, 'Position', get(0, 'Screensize'));
end
hold on;


%% ---------------- Data --------------- %%
chrData = obj.chrFDictionary;
CNReference  = obj.CNReference;
binSize = obj.binSize;
%
chrs = obj.targetChrs;
if (chrs == 0) 
    chromosomes = [1:1:length(obj.chrNames)];
else
    chromosomes = chrs;
end
noChrs = length(chromosomes);

%% ------------ Processing --------------%%
genomeRD = [];
%%%
for i=1:noChrs
    targetChrIndex = chromosomes(i);
    %%%%%%%%%% Data Reading %%%%%%%%

    orgData = chrData(targetChrIndex);
    copyNumber = CNReference;
    normData = orgData*2/copyNumber;
    %
    clippedDataLen = length(normData);
    clippedData    = normData;

    %%%%%%%%%% Scaling %%%%%%%%%
    n = round(nBinSize/binSize);
    %
    a = clippedData;
    sClippedData = arrayfun(@(k) mean(a(k:k+n-1)),1:n:length(a)-n+1)';
    %
    genomeRD = [genomeRD; sClippedData];
end
%%%
genomeRD = genomeRD(genomeRD <8);
h = histogram(genomeRD,200);
plot([2 2],[0 max(h.Values)],'k', 'LineWidth',2);
%
ylabel('Frequency', 'FontSize', 8,'FontWeight','bold');
xlabel('Copy number', 'FontSize', 8,'FontWeight','bold');
title('Genome histogram', 'FontSize', 10,'FontWeight','bold');
hold off;


%%---------------------- Save figure ---------------------------%%
dir = obj.outputDirectory;
if(exist(dir,'dir') ~= 7)
    mkdir(dir);
end
[~,objFileName,~] = fileparts(obj.bamFile);
%
if(saveResult == 1)
    ff = strcat('-f102');
    hh = strcat(obj.outputDirectory, '/', objFileName, '_GenomeWide_histogram');
    %
    print(ff,hh,'-dpng');
    savefig(hh);
    close All;
end


end
%%%