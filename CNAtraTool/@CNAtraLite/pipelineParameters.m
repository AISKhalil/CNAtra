function obj = pipelineParameters(obj)
%Estimating the CNAtra thresholds based on data coverage using exponential regression model.

dataCoverage = obj.dataCoverage;
%%%%----- Regression Model (Analysis parameters vs data coverage) -----%%
%++1) resolution: minimum alteration region that can be detected accurately based on the data coverage. It is also used as Savitzky frame length.
%%Power-decay model
a = 38.99;
b = -0.3772;
c = 17.08; 
d = -0.008899;
x = dataCoverage;
minRegionWidthForTestInBins = 15;
obj.resolution = max(round(a*exp(b*x) + c*exp(d*x)),minRegionWidthForTestInBins);


%++2) coverageThreshold: threshold for considering alteration region as practically significant region
%% Amplification
a2 = 0.3434;
b2 = -0.8585;
c2 = 0.6872;
d2 = -0.02144;
x2 = dataCoverage;
obj.amplificationThreshold = min(max(round(a2*exp(b2*x2) + c2*exp(d2*x2),3), 0.5),1);

%% Deletion
a3 = 0.5295;
b3 = -1.14;
c3 = 0.6581;
d3 = -0.01981;
x3 = dataCoverage;
obj.deletionThreshold = min(max(round(a3*exp(b3*x3) + c3*exp(d3*x3),3), 0.5),1);

