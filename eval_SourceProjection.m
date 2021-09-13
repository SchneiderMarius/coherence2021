% Model of Figure S3
% add fieldtrip toolbox to your path
clear all
addpath(fullfile(cd,'functions'));

%%
ModulationStrength = logspace(-2,0,20)*0.5;

frate   = cell(1,length(ModulationStrength));
modstr  = cell(1,length(ModulationStrength));
Ntest   = cell(1,length(ModulationStrength));
stepnum = cell(1,length(ModulationStrength));
plt     = cell(1,length(ModulationStrength));

for cnt1 = 1 : length(ModulationStrength)
    frate{cnt1}     = 2;
    modstr{cnt1}    = ModulationStrength(cnt1);
    Ntest{cnt1}     = 10000;
    stepnum{cnt1}   = 20;
    plt{cnt1}       = 0;
end

[result] = cellfun(@ProjectionSource,frate,modstr,Ntest,stepnum,plt);    

save('Fig3C','result')