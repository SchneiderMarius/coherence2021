% Model of Figure S6
% add fieldtrip toolbox to your path
clear all
addpath(fullfile(cd,'functions'));
load(fullfile(cd,'par','ARparameter'));

%%

par             = parAR;
par.fsample     = 1000;
par.numTrial    = 800;
par.time        = 4;    % Time in sec
par.N           = 500;  % Number of Neurons
par.cw          = 0.05; % projection strength
par.fac         = 0.5;  % factor determining SOS
par.g           = 1;    % gamma parameter from SSM model

result = LIFmod(par);
