%% Model of Figure S7
% add fieldtrip toolbox to your path
clear all
addpath(fullfile(cd,'functions'));

%% 
load(fullfile('par','Poissonparameter'));

numTrial        = 500;
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
            result.snr(cnt2,cnt3) = out(id).snr_temp;
            result.CohSim(cnt2,cnt3) = out(id).Cohf;  
            result.cohAna(cnt2,cnt3) = sqrt(1./(1+(1./(cw(cnt3).^2*(1+out(id).snr_temp)))));
            id = id + 1;
        end
    end
    
    save(fullfile(cd,sprintf('LFPpart_%s',func{cnt1})),'result');
    clear result
end

%% 2) Spike Coherence
clear all
addpath(fullfile(cd,'functions'));

Test.maxfire    = 2;
Test.cw         = [logspace(-3,0,15)];
Test.alpha      = [logspace(-2,0,10) 2 4 6 8 10 14 50 100];
Test.g          = 1;
Test.fac        = fliplr(Test.alpha);
Test.N          = 1000;
steps           = 20;

id = 1;        
for cnt1 = 1 : length(Test.fac)
    for cnt2 = 1 : length(Test.cw)
        cw{id} = Test.cw(cnt2);
        alpha{id} = Test.fac(cnt1);
        N{id} = Test.N;
        maxfire{id} = Test.maxfire;
        g{id} = Test.g(1);
        stps{id} = steps; 
        id = id + 1;
    end
end
out = cellfun(@PoissonSpikeCoh,cw,alpha,maxfire,g,N,stps);      

id = 1;        
for cnt3 = 1 : length(Test.fac)
    for cnt4 = 1 : length(Test.cw)
        result.Frate(cnt3,cnt4) = out(id).Frate;
        result.cw(cnt3,cnt4) = Test.cw(cnt4);
        result.snr(cnt3,cnt4) = out(id).SNR;               
        for cnt5 = 1 : steps
            result.pow(cnt3,cnt4,cnt5,:,:) = out(id).pow(:,:,cnt5);
            result.CohSim(cnt3,cnt4,cnt5) = out(id).Cohf(cnt5);  
        end
        id = id + 1;                
    end
end
result.strcoh   = 'snr_cw_N';
result.strpow   = 'snr_cw_N_Signal_frq';
result.maxfire  = Test.maxfire(cnt2);
result.N        = out(1).N;
result.powfrq   = out(1).f;

save(fullfile(cd,'SpikeCoh'),'result')
clear result

%% 3) Spike PPC
clear all
addpath(fullfile(cd,'functions'));

Test.cw         = logspace(-3,0,15);
Test.alpha      = [logspace(-2,0,10) 2 4 6 8 10 14 50 100];
Test.fac        = fliplr(Test.alpha);
Test.g          = 1;
Test.maxfire    = 1;
Test.N          = 50;
steps           = 1;
Test.cw         = Test.cw(9);
Test.fac        = Test.fac(9);

id = 1;        
for cnt1 = 1 : length(Test.fac)
    for cnt2 = 1 : length(Test.cw)
        cw{id}      = Test.cw(cnt2);
        alpha{id}   = Test.fac(cnt1);
        N{id}       = Test.N;
        maxfire{id} = Test.maxfire;
        g{id}       = Test.g(1);
        stp{id}     = steps(1);
        ind{id}     = id; 
        
        id = id + 1;
    end
end

out = cellfun(@PoissonSpikePPC,cw,alpha,maxfire,g,N,stp);

id = 1;   

result.PPCfreq = out(1).PPCfrq;
result.Cohfreq = out(1).Cohfrq;
for cnt1 = 1 : length(Test.fac)
    for cnt2 = 1 : length(Test.cw)
        result.Frate(cnt1,cnt2) = out(id).Frate;
        result.cw(cnt1,cnt2) = Test.cw(cnt2);
        result.snr(cnt1,cnt2) = out(id).SNR;
        result.PPC{cnt1,cnt2} = out(id).PPC;  
        result.Coh{cnt1,cnt2} = out(id).Coh;  
        id = id + 1;                
    end
end
result.str      = 'snr_cw';
result.maxfire  = Test.maxfire;
result.N        = Test.N;
save(fullfile('SpikePPC'),'result')
