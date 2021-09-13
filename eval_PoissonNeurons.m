%% Model of Figure S7
% add fieldtrip toolbox to your path
clear all
addpath(fullfile(cd,'functions'));

%% 
load(fullfile('par','Poissonparameter'));
numTrial        = 2%500;
num             = 2;
params.fsample = 1000;
params.rate     = 1/params.fsample;
params.time     = 10;
params.delay    = 4;
params.a        = 2/3;
params.freq     = [17 19];
len             = 1/params.rate*params.time;
time            = params.rate:params.rate:params.time;
cw              = [logspace(-3,0,15)];
SNR             = [logspace(-2,0,10) 2 4 6 8 10 14 50 100];
fac             = fliplr(SNR);
func            = {'sigmoid','linear'};

numSim = length(SNR)*length(cw);

%
for cnt1 = 1 : length(func)
    result.func = func{cnt1};
    
    id = 1;
    funcin = cell(1,numSim);
    parARin = cell(1,numSim);
    paramsin = cell(1,numSim);
    numTrialin = cell(1,numSim);
    facin = cell(1,numSim);
    cwin = cell(1,numSim);
    for cnt2 = 1 : length(SNR)
        for cnt3 = 1 : length(cw)
            funcin{id} = func{cnt1};
            parARin{id} = parAR;
            paramsin{id} = params;
            numTrialin{id} = numTrial;
            facin{id} = fac(cnt2);
            cwin{id}  = cw(cnt3);
            id = id + 1;
        end
    end
    out = cellfun(@SSM2,parARin,paramsin,numTrialin,facin,cwin,funcin);                

    id = 1;
    for cnt2 = 1 : length(SNR)
        for cnt3 = 1 : length(cw)
            result.cw(cnt2,cnt3) = cw(cnt3);
            result.snr(cnt2,cnt3) = out{id}.snr_temp;
            result.CohSim(cnt2,cnt3) = out{id}.Cohf;  
            result.cohAna(cnt2,cnt3) = sqrt(1./(1+(1./(cw(cnt3).^2*(1+out{id}.snr_temp)))));
            id = id + 1;
        end
    end
    
    save(fullfile(cd,sprintf('LFPpart_%s',func{cnt1})),'result');
    clear result
end

%% 2) Spike Coherence
clear all
addpath /mnt/hpx/opt/ESIsoftware/matlab/slurmfun
addpath /opt/ESIsoftware/matlab/slurmfun
addpath(fullfile(cd,'helper_functions_model'));

%save_folder = fullfile(cd,'results','LFPpartModCoh4');
SaveDir = fullfile(cd,'results');

% original:
%Test.func = {'sigmoid1','sigmoid2','sigmoid3','linear'};
%Test.maxfire = [1 2 10];%round([logspace(0,2,10)]);%[2 round([ logspace(1,4,10)]/10)*10];
% plotted
%nm = {'Modpart_sigmoid2_rate2.mat','Modpart_sigmoid3_rate2.mat'}
%Test.func = {'sigmoid2','sigmoid3'};
%Test.maxfire = [2];%round([logspace(0,2,10)]);%[2 round([ logspace(1,4,10)]/10)*10];
Test.maxfire    = 2;
Test.cw         = [logspace(-3,0,15)];
Test.alpha      = [logspace(-2,0,10) 2 4 6 8 10 14 50 100];%round([logspace(-1,0,10) 2 4]*100)/100;% 14 50 100];
Test.g          = 1;%[0 0.5 0.95 1];
Test.fac        = fliplr(Test.alpha);
Test.N          = 1000;
steps           = 20;

%for cnt2 = 1 : length(Test.maxfire)
id = 1;        
for cnt3 = 1 : length(Test.fac)
    for cnt4 = 1 : length(Test.cw)
        cw{id} = Test.cw(cnt4);
        alpha{id} = Test.fac(cnt3);
        N{id} = Test.N;
        maxfire{id} = Test.maxfire;
        g{id} = Test.g(1);
        fl{id} = SaveDir;
        stps{id} = steps; 
        id = id + 1;
    end
end
out = slurmfun(@PoissonSpikeCoh,cw,alpha,maxfire,g,N,stps,fl,'partition', '8GBXS');      
%        out = cellfun(@SpikeModPartTransferCohtestRate,cw,alpha,maxfire,g,N,stps,func)

id = 1;        
for cnt3 = 1 : length(Test.fac)
    for cnt4 = 1 : length(Test.cw)
        result.Frate(cnt3,cnt4) = out{id}.Frate;
        result.cw(cnt3,cnt4) = Test.cw(cnt4);
        result.snr(cnt3,cnt4) = out{id}.SNR;               
        for cnt5 = 1 : steps
            result.pow(cnt3,cnt4,cnt5,:,:) = out{id}.pow(:,:,cnt5);
            result.CohSim(cnt3,cnt4,cnt5) = out{id}.Cohf(cnt5);  
        end
        id = id + 1;                
    end
end
result.strcoh = 'snr_cw_N';
result.strpow = 'snr_cw_N_Signal_frq';

result.maxfire = Test.maxfire(cnt2);
result.N = out{1}.N;
result.powfrq = out{1}.f;

save(fullfile(SaveDir,'FigS7','SpikeCoh'),'result')
%save(fullfile(save_folder,sprintf('Modpart_%s_rate%d',Test.func{cnt1},result.maxfire)),'result')
clear result
%end



%% 3) Spike PPC
%% b) PPC
clear all
addpath /mnt/hpx/opt/ESIsoftware/matlab/slurmfun
addpath /opt/ESIsoftware/matlab/slurmfun
addpath(fullfile(cd,'helper_functions_model'));

save_folder = fullfile(cd,'results','FigS7');
SaveDir = fullfile(cd,'results');
%save_folder = fullfile(cd,'results','LFPpartModPPC4');
% original:
%Test.func = {'sigmoid1','sigmoid2','sigmoid3','linear'};
%Test.maxfire = [1 2 10];%round([logspace(0,2,10)]);%[2 round([ logspace(1,4,10)]/10)*10];

% plotted
%nm = {'Modpart_sigmoid2_rate2.mat','Modpart_sigmoid3_rate2.mat'}
%Test.func = {'sigmoid2'};
%Test.maxfire = [1,2];%round([logspace(0,2,10)]);%[2 round([ logspace(1,4,10)]/10)*10];



Test.cw = [logspace(-3,0,15)];
Test.alpha = [logspace(-2,0,10) 2 4 6 8 10 14 50 100];%round([logspace(-1,0,10) 2 4]*100)/100;% 14 50 100];
% Test.cw = 1;
% Test.alpha = [logspace(-2,0,10) 2];%round([logspace(-1,0,10) 2 4]*100)/100;% 14 50 100];
%Test.maxfire = [1 2 10];%round([logspace(0,2,10)]);%[2 round([ logspace(1,4,10)]/10)*10];
%Test.g = 1;%[0 0.5 0.95 1];
%Test.func = {'sigmoid','sigmoid2','sigmoid3','linear'};

Test.fac= fliplr(Test.alpha);
Test.g = 1;%[0 0.5 0.95 1];
Test.maxfire = 1;
Test.N       = 50;
steps        = 1;
Test.cw      = Test.cw(9);
Test.fac     = Test.fac(9)

%for cnt1 = 1 : length(Test.func)
%    for cnt11 = 1 : length(Test.N)
%        for cnt2 = 1 : length(Test.maxfire)
            %result.func = Test.func{cnt1};
id = 1;        
for cnt3 = 1 : length(Test.fac)
    for cnt4 = 1 : length(Test.cw)
        cw{id} = Test.cw(cnt4);
        alpha{id} = Test.fac(cnt3);
        N{id} = Test.N;
        maxfire{id} = Test.maxfire;%(cnt2);
        g{id} = Test.g(1);
        folder{id} = SaveDir;
        stp{id} = steps(1);
        ind{id} = id; 
        id = id + 1;
    end
end

out = slurmfun(@PoissonSpikePPC,cw,alpha,maxfire,g,N,stp,folder,'partition', '8GBXS')

% if Test.maxfire(cnt2)>5
% %                 out = slurmfun(@SpikeModPartTransferPPC,cw,alpha,maxfire,g,N,stp,func,'partition', '16GBS')
%     %PPC from SuperNeuron
%     out = slurmfun(@PoissonSpikePPC,cw,alpha,maxfire,g,N,stp,folder,'partition', '16GBS')            
% else
% %                 out = slurmfun(@SpikeModPartTransferPPC,cw,alpha,maxfire,g,N,stp,func,'partition', '8GBXS')
%     %PPC from SuperNeuron    
%     out = slurmfun(@PoissonSpikePPC,cw,alpha,maxfire,g,N,stp,folder,'partition', '8GBXS')
% end
%             out = cellfun(@SpikeModPartTransferPPC2,cw,alpha,maxfire,g,N,ind,func)

id = 1;    
result.PPCfreq = out{1}.PPCfrq;
result.Cohfreq = out{1}.Cohfrq;

for cnt3 = 1 : length(Test.fac)
    for cnt4 = 1 : length(Test.cw)
        result.Frate(cnt3,cnt4) = out{id}.Frate;
        result.cw(cnt3,cnt4) = Test.cw(cnt4);
        result.snr(cnt3,cnt4) = out{id}.SNR;
        result.PPC{cnt3,cnt4} = out{id}.PPC;  
        result.Coh{cnt3,cnt4} = out{id}.Coh;  
        id = id + 1;                
    end
end
result.str = 'snr_cw';
result.maxfire = Test.maxfire;%(cnt2);
%result.func = Test.func{cnt1};
result.N = Test.N;%(cnt11);

save(fullfile(SaveDir,'FigS7','SpikePPC'),'result')
%save(fullfile(save_folder,sprintf('ModpartPPC_%s_rate%d',Test.func{cnt1},cnt2)),'result')
clear result
%        end
%    end
%end






