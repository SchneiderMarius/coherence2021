% Model of Figure S7
% add fieldtrip toolbox to your path

clear all
addpath(fullfile(cd,'functions'));

load(fullfile(cd,'par','ARparameter'));

par = parAR;
par.numTrial    = 500;
par.N           = 1000; % number of neurons
par.fsample     = 1000; 
par.time        = 10;   % Time in sec
par.stepnum     = 20;   % evaluated population size 
par.frate       = 2;    % firing rate    
par.g           = 1;    % parameter gamma from SSM2
par.a           = 2/3;  % slope of 1/f
par.cw          = 0.05; % projection strength
par.fac         = 5;    % determines SOS
par.modstr      = 0.1;  % determines Phase-Locking

result = PoissonSpikeCoh(par);      
