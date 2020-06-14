classdef CNAtraLite < handle
    properties (GetAccess=public, SetAccess=public)

		%%%% -------------- Input-files ----------------%%
		
		%bamFile: Input bam-formatted file to be analyzed 
		bamFile

		%outputDirectory: Output directory to save results and figures at.
		outputDirectory = './CNAtraResults';



		%%%% ------------ RD-calculator  -------------%%

		%memoryFootPrint: Number of reads Per iterations for GC-calculations from the input reads. User can choose it based on RAM size.
		memoryFootPrint = 100000;

		%MAPQ_Score: A threshold for filtering the reads from the BAM file. Default value is 1 assuming Bowtie2 mapping, however user can change it based on his short sequence mapper.
		MAPQ_Score = 1;


		%%%% -------------- CNV-caller --------------%%
		
		%ploidyLevel: Expected number of chromosome sets (default = free), user can set it based on input-cell ploidy {'free', 'diploid', 'triploid', 'tetraploid'}.
		ploidyLevel = 'free'; 
        		
		%minimumIBsize: A minimum segment size to be considered as Iso-copy numeric block (in Bins), default size = 1Mb.
		minimumIBsize = 1000;

        %resolution: Minimum bin size that is computed using regression model based on the data coverage. We used it as a frame length for Savitzky Golay filter.
		resolution;

		%amplificationThreshold: Threshold for identifying the class2 amplifications.
		amplificationThreshold;

		%deletionThreshold: Threshold for identifying the class2 deletions.
		deletionThreshold;
		
		%minAlterationRank: Minimum rank of an alteration to be considered as significant FA: (2) Class1 & Class2, (5) Class2 only.
	    minAlterationRank = 5;
		

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
        
		%chrNames: chromosome names.
		chrNames;
		
		%chrLengths: chromosome lengths.
		chrLengths;

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

		%removeShortRegions: Filter short regions based on the data resolution. 
		removeShortRegions = 1;


	end
	%%%


	%%%%%%%
	%%%%%%%
	methods
        %%%% ------------ Constructor ------------%%%%%
        function obj = CNAtraLite(rdFile)
        %A Read-Depth Based Hierarchical Approach to Detect Copy Number Variations in Low-Coverage and Complex Copy Number Profiles.
	      if nargin == 1
	         obj.bamFile = rdFile;
	      else
	         error('No enough inputs')
	      end
        end    
		

        %%%% ------- RD-calculator ------- %%%%

        RDcalculator(obj)
        bamBinning(obj, rdFile)
        pipelineParameters(obj)
		

        %%%% --------- CNV-caller -------- %%%%

		CNVcaller(obj)
		[CNReference] = copyNumberReference(obj, signal)
		[segmentStart, segmentEnd, largeSegmentStart, largeSegmentEnd] = rdSavitzkySegmentation(obj, signal)
		[finalSegmentStart, finalSegmentEnd] = shortSegmentsMerging(obj, orgSignal, orgSegmentStart, orgSegmentEnd, CNReference)
		[segmentsInfo, regionsPerSegment] = CNRegionCalling(obj, orgSignal, orgSegmentStart, orgSegmentEnd, finalSegmentStart, finalSegmentEnd, CNReference)
		[finalSegmentInfo,impRegionsPerSegment] = CNRegionFiltering(obj, signal, targetChrIndex, segmentsInfo, regionsPerSegment)
		[segmentsData, regionsData] = CNResultsUpdate(obj, targetChrIndex, segmentsInfo, impRegionsPerSegment)


		%%%% ------ Visualization ------- %%%%
		%Plot chromosome or a specific area
		CNAPlot(obj, saveResult, varargin)
		CNVsTrackPlot(obj, saveResult, varargin) 
		
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
