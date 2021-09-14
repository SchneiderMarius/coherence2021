% Model of Figure S3
% add fieldtrip toolbox to your path
clear all
addpath(fullfile(cd,'functions'));
load(fullfile(cd,'par','ARparameter'));

%%

ModulationStrength  = logspace(-2,0,20)*0.5;
par                 = parAR;
par.numTrial        = 800;
par.N               = 100;
par.fsample         = 1000;
par.time            = 3;
par.frate           = 2;
par.stepnum         = 20; 

parameter = cell(1,length(ModulationStrength));

for cnt1 = 1 : length(ModulationStrength)
    parameter{cnt1}         = par;
    parameter{cnt1}.modstr	= ModulationStrength(cnt1);
end

result = cellfun(@ProjectionSource,parameter);    
