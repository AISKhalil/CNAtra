classdef CNAtra < handle
    properties (GetAccess=public, SetAccess=public)

		%%%% -------------- Input-files ----------------%%
		
		%bamFile: Input bam-formatted file to be analyzed 
		bamFile

		%gcWindsMatFile: File that contains the pre-calculated gc-contents per bin. A default hg19 file "CNAClassDirectory+/annotations/hg19/AnshulMappability/gcWinds.readLength100.mat" 
		gcWindsMatFile;

		%mapMatFile: File that includes mappability scores of the genome (at 1Kb bins) which can be to filter highly repeated or unmappable regions. A default hg19-based file "CNAClassDirectory+/annotations/hg19/ChrisMillerGCs/mapTracks.hg19.101.mat" is included by default if user doesn't specify one. This file can be used for Single-end reads with read length >= 100bps or Paired-end reads. For other reads, it is better to load other files (see tool documentations)	
		mapMatFile;

		%gapFile: File that includes unmappable hg19 regions. A default hg19-based file "CNAClassDirectory+/annotations/hg19/hg19gaps.bed" is included by default if user doesn't specify one.		
		gapFile;

		%blackListFile: File that includes black-listed regions. A default hg19-based file "CNAClassDirectory+/annotations/hg19/hg19blackListed.bed" is included by default if user doesn't specify one.
		blackListFile;	

		%centromeresFile: File that includes locations of centromeres. A default hg19 file "CNAClassDirectory+/annotations/hg19/hg19centromeres.bed" is included by default if user doesn't specify one.
		centromeresFile;

		%telomeresFile: File that includes locations of telomeres. A default hg19 file "CNAClassDirectory+/annotations/hg19/h19telomeres.bed" is included by default if user doesn't specify one.
		telomeresFile;

		%outputDirectory: Output directory to save results and figures at.
		outputDirectory = './CNAtraResults';



		%%%% ------------ RD-calculator  -------------%%

		%memoryFootPrint: Number of reads Per iterations for GC-calculations from the input reads. User can choose it based on RAM size.
		memoryFootPrint = 100000;

		%MAPQ_Score: A threshold for filtering the reads from the BAM file. Default value is 1 assuming Bowtie2 mapping, however user can change it based on his short sequence mapper.
		MAPQ_Score = 1;

		%gcCorrectionMethod: Method that is used to correct the read-depth value per bin. 0:No GC-correction is applied, 1: Genome-wise GC-correction using pre-calculated ChrisMiller gcTracks, 2: Genome-wise GC-correction using gcTracks that are calculated from the input data.
        gcCorrectionMethod = 1;



		%%%% -------------- CNV-caller --------------%%
		
		%ploidyLevel: Expected number of chromosome sets (default = free), user can set it based on input-cell ploidy {'free', 'diploid', 'triploid', 'tetraploid'}.
		ploidyLevel = 'free'; 
        		
		%minimumIBsize: minimum segment size to be considered as Iso-copy numeric block (in Bins), default size = 1Mb.
		minimumIBsize = 1000;

        %resolution: Minimum bin size that is computed using regression model based on the data coverage. We used it as a frame length for Savitzky Golay filter.
		resolution;

		%amplificationThreshold: Threshold for identifying the class2 amplifications.
		amplificationThreshold;

		%deletionThreshold: Threshold for identifying the class2 deletions.
		deletionThreshold;
		
		%minAlterationRank: Minimum rank of an alteration to be considered as significant FA: (2) Class1 & Class2, (5) Class2 only.
	    minAlterationRank = 2;
		
		%minMappabilityThreshold: Minimum mappability-score of a bin to be kept for CNAtra analysis.(if 1: minMappabilityThreshold is computed as 10thPercentile of the mappability scores based on read length, else: the minMappabilityThreshold is used as a threshold).
		minMappabilityThreshold = 0.5;

		%maximumFalseBinsAllowed: Maximum allowed region percentage to keep alteration regions
		maximumFalseBinsAllowed = 0.5;
        


		%% ------------ Data-statistics ------------%%
		
		%numberOfReads: Total number of reads in the input file;
		numberOfReads;

		%dataCoverage: Genome coverage of the dataset = number_of_reads * read_length/genome_length.
	    dataCoverage;
		
		%readLength: Read length in terms of bps.
		readLength;

		%readType: Single-end/paired-end
		readType; 



		%% ------------ Data Variables ------------%%

		%chrRawDictionary: ReadDepth signal per chromosome before GC-correction, filtering black-listed, centromeres, and telomeres.
		chrRawDictionary;

		%chrDictionary: ReadDepth signal per chromosome after GC-correction.
		chrDictionary;

		%chrFDictionary: ReadDepth signal per chromosome after filtering black-listed, centromeres, and telomeres.
		chrFDictionary;

		%chrFIndex: Dictionary of filtered bins;
		chrFIndex;
        
		%chrNames: chromosome names.
		chrNames;
		
		%chrLengths: chromosome lengths.
		chrLengths;

		%chrMappabilityTracks: Dictionary of mappability tracks per chromosome.
		chrMappabilityTracks;

		%chrGCTracks: Dictionary of GC tracks per chromosome.
		chrGCTracks;

		%chrCentroTeloBoundaries: Dictionary of extended telomeres and centromeres per chromosome.
		chrCentroTeloBoundaries;

		%chrBlackListedBoundaries: Dictionary of extended black-listed per chromosome.
		chrBlackListedBoundaries;

		%chrGapBoundaries: Dictionary of gap regions per chromosome.
		chrGapBoundaries;

		%CNReference: Copy number reference.
		CNReference;

		%segmentsInfoDic: Dictionary of IBs per chromosome.
		segmentsInfoDic;

		%regionsInfoDic: Dictionary of CNVs per chromosome.
		regionsInfoDic;	

    end
	%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	properties (GetAccess=private, SetAccess=private)

		%targetChrs: Indices of the given-chromosomes: chr1:1, chr2:2, ....., chrX:23.
        targetChrs;

		%binSize: Number of base pairs per bin. 
        binSize = 1000;

		%sampleSize: Sample-size that is used for computing the p-value of the t-test (default value = 35).
		sampleSize = 35;

		%noIterations: Number of iterations to repeat sampling t-test (default value = 10).
        noIterations = 10;

		%alpha: Threshold for t-test to identify the statistically significant segments (default value = 0.05).
        alpha = 0.05;

		%keepCNVinBCT: Flag to keep regions around black-listed regions. %0: remove any region around black-listed, centromeres, and telomeres.
	    keepCNVinBCT = 0;      

		%removeLowMappabilityBins: Flag to remove low-mappability bins.		
		removeLowMappabilityBins = 1;

		%removeShortRegions: Filter short regions based on the data resolution. 
		removeShortRegions = 1;

		%uncertaintyDistanceToCentro: Number of bins, surrounding centromeres with normal RD-value to extend the centromere boundaries (default size: 100Kb). 
		uncertaintyDistanceToCentro = 100;
		
		%uncertaintyDistanceToTelo: Number of bins, surrounding telomeres with normal RD-value to extend the telomeres boundaries (default size: 10Kb). 
		uncertaintyDistanceToTelo = 10;
		
		%chrInputGCContents: Dictionary of GC tracks per chromosome.
		chrInputGCContents;

		%intialSegmentsInfoDic:Dictionary of initial segments.
		intialSegmentsInfoDic;

	end
	%%%


	%%%%%%%
	%%%%%%%
	methods
        %%%% ------------ Constructor ------------%%%%%
        function obj = CNAtra(rdFile,CNAClassDirectory)
        %A Read-Depth Based Hierarchical Approach to Detect Copy Number Variations in Low-Coverage and Complex Copy Number Profiles.
	      if nargin == 2
	         obj.bamFile = rdFile;
	         obj.blackListFile    = strcat(CNAClassDirectory,'/annotations/hg19/hg19blackListed.bed');
	         obj.centromeresFile  = strcat(CNAClassDirectory,'/annotations/hg19/hg19centromeres.bed');
	         obj.telomeresFile    = strcat(CNAClassDirectory,'/annotations/hg19/hg19telomeres.bed');
	         obj.gapFile          = strcat(CNAClassDirectory,'/annotations/hg19/hg19gaps.bed');
	         obj.mapMatFile       = strcat(CNAClassDirectory,'/annotations/hg19/AnshulMappability/mapTracks.hg19.101.mat');
	         obj.gcWindsMatFile   = strcat(CNAClassDirectory,'/annotations/hg19/ChrisMillerGCs/gcWinds.hg19.readLength100.mat');
	      else
	         error('No enough inputs')
	      end
        end    
		

        %%%% ------- RD-calculator ------- %%%%

        RDcalculator(obj)
        bamBinning(obj, rdFile)
        RDcorrection(obj)
        binFiltering(obj)
        pipelineParameters(obj)
		

        %%%% --------- CNV-caller -------- %%%%

		CNVcaller(obj)
		[CNReference] = copyNumberReference(obj, signal)
		[segmentStart, segmentEnd, largeSegmentStart, largeSegmentEnd] = rdSavitzkySegmentation(obj, signal)
		[finalSegmentStart, finalSegmentEnd] = shortSegmentsMerging(obj, orgSignal, orgSegmentStart, orgSegmentEnd, CNReference)
		[segmentsInfo, regionsPerSegment] = CNRegionCalling(obj, orgSignal, orgSegmentStart, orgSegmentEnd, finalSegmentStart, finalSegmentEnd, CNReference)
		[finalSegmentInfo,impRegionsPerSegment] = CNRegionFiltering(obj, signal, targetChrIndex, segmentsInfo, regionsPerSegment, refIndices, centroTeloBoundaries)
		[segmentsData, regionsData] = CNResultsUpdate(obj, targetChrIndex, segmentsInfo, impRegionsPerSegment, centroTeloBoundaries)
		[segmentsInfo] = neighborSegmentsMerging(obj, targetChrIndex, orgSignal, orgSegmentStart, orgSegmentEnd, CNReference)	


		%%%% ------ Visualization ------- %%%%

		%Plot chromosome or a specific area
		CNAPlot(obj, saveResult, varargin)

		%Plot amplification or deletion region by its number
		CNARegionPlot(obj, saveResult, varargin)

		%Plot CNV track
		CNVsTrackPlot(obj, saveResult, varargin) 

		%Plot Initial-segments
		CNAPlotInitialSegments(obj, saveResult, varargin)

		%Estimate the copy number of a region 
		[regionCN] = CNAEstimator(obj, varargin)

		%Ploidy-Test
		ploidyTest(obj, saveResult)

		%Plot genome RD signal
		plotGenome(obj, saveResult, nBinSize)

		%Plot genome RD histogram
		plotGenomeHistogram(obj, saveResult, nBinSize)

    end
    %%% 

end   
%%%