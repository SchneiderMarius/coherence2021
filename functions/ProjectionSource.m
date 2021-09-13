function  [out] = ProjectionSource(frate,modstr,Ntest,stepnum,plt,folder)

if nargin<1
    frate = 2;
end
if nargin<2
   modstr = 0.1; 
end
if nargin<3 || isempty(Ntest)
   Ntest = 10000; 
end
if nargin<4 || isempty(stepnum)
    stepnum = 20;
end
if nargin<5 || isempty(plt)
    plt = 0;
end


load(fullfile(cd,'par','ARparameter'));

%% noise
Nsub            = 100;
numTrial        = 400;
params.fsample  = 1000;
num             = 2;
params.rate     = 1/params.fsample;
params.time     = 3;
len             = 1/params.rate*params.time;
time            = params.rate:params.rate:params.time;

for cnt0 = 1 : numTrial
    S{cnt0} = zeros(length(time),1);
    for cnt1 = 1 : size(parAR{1}.coeff,1)
        S{cnt0}(1:2,cnt1) =  randn(2,1);
        for cnt2 = 3 : length(time)              
           S{cnt0}(cnt2,cnt1) = parAR{1}.coeff(cnt1,:) * S{cnt0}(cnt2-2:cnt2-1,cnt1) + randn;         
        end
    end
end

params.Ntest = Ntest;

if plt == 1
    stp = [Ntest/100, Ntest];
    out.Signal.Mod = S{1};
    [spikeMat, Frate,LF{1}]           = inhomopp(S{1},modstr,params,frate,stp); 
    
    out.Signal.spikeMat = spikeMat(1:1000,:);
    clear spikeMat
    LF{1} = full(LF{1});
    
    [~,td,~] = fieldtrip_it(LF,params);
    cfg          = [];
    cfg.lpfilter = 'yes';
    cfg.lpfreq   = 100;
    cfg.hpfilter = 'yes';
    cfg.hpfreq   = 1;
    [td]         = ft_preprocessing(cfg,td);    
    
    out.Signal.Source     = td.trial{1}(2,:);
    out.Signal.Projection = td.trial{1}(1,:);
    out.Signal.time       = td.time{1}; 
end

LF = cell(1,length(S));

SP = cell(1,length(S));
for cnt1 = 1 : length(S)  
    [SPmat, Frate(cnt1),LF{cnt1}]     = inhomopp(S{cnt1},modstr,params,frate,steps);        
    SP{cnt1} = SPmat(1:Nsub,:);
end

LFP = cell(1,length(S));
for cnt1 =1 : length(S)
    LFP{cnt1} = [S{cnt1}';full(SP{cnt1})]; 
    SP{cnt1} = [];
end
[~,tdSP,~] = fieldtripIT(LFP,params);

% evalulate PPC
cfg              = [];
cfg.method       = 'mtmconvol';
cfg.foi          = [1:1:100];
cfg.timwin       = [-0.4 0.4];
cfg.taper        = 'hanning';
cfg.t_ftimwin    = ones(length(cfg.foi),1);
cfg.spikechannel = tdSP.label(2:Nsub+1);
cfg.channel      = tdSP.label(1);
cfg.toi          = 0:0.1:tdSP.time{1}(end);
stsFFT           = ft_spiketriggeredspectrum(cfg, tdSP);    

PPC_temp = [];
for cnt2 = 1 : Nsub
    cfg               = [];
    cfg.method        = 'ppc0';
    cfg.spikechannel  = tdSP.label{cnt2+1};
    cfg.channel       = tdSP.label{1};
    cfg.avgoverchan   = 'unweighted';
    statSts1          = ft_spiketriggeredspectrum_stat(cfg,stsFFT);   
    PPC_temp(cnt2,:)  = statSts1.ppc0;     
end

out.PPCfrq           = statSts1.freq;     
out.PPC(1,:)         = mean(PPC_temp,1);    

% evalulate Coh
cfg              = [];
cfg.output       = 'fourier';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          =  [1:100];
cfg.t_ftimwin    = ones(length(cfg.foi),1); 
cfg.toi          = 0:0.1:tdSP.time{1}(end);
ff               = ft_freqanalysis(cfg, tdSP);    

chlcmb = {};
lbl = unique(tdSP.label);
chlcmb(1:length(lbl)-1,1) = lbl(1);
chlcmb(1:length(lbl)-1,2) = lbl(2:end);

coh_temp = [];

for cnt1 = 1 : length(chlcmb)
    cfg             = [];
    cfg.channelcmb  = chlcmb(cnt1,:);
    cfg.method      = 'coh';
    coh             = ft_connectivityanalysis(cfg, ff); 

    cfg              = [];
    cfg.avgovertime  = 'yes';
    cfg.nanmean      = 'yes';
    coh              = ft_selectdata(cfg,coh); 

    coh_temp = cat(1,coh_temp,coh.cohspctrm);
end

out.SPLFP(1,:)    = mean(coh_temp,1);
out.SPLFPfrq      = coh.freq;
out.Ntotal        = Ntest;
out.Nproj         = steps;
out.g             = steps/Ntest;

[val,~]  = max(squeeze(out.SPLFP),[],2); 
out.SPLFPPeak     = val;
%% 
for cnt1 = 1 : length(steps)-1
    for cntT =1 : length(LF)
        LFP{cntT} = [full(LF{cntT}(cnt1,:));full(LF{cntT}(end,:))]; 
    end
    
    [~,td,~] = fieldtrip(LFP,params);
    
    cfg          = [];
    cfg.lpfilter = 'yes';
    cfg.lpfreq   = 100;
    cfg.hpfilter = 'yes';
    cfg.hpfreq   = 1;
    [td]         = ft_preprocessing(cfg,td);
    
    cfg              = [];
    cfg.output       = 'fourier';
    cfg.method       = 'mtmconvol';
    cfg.taper        = 'hanning';
    cfg.foi          =  [1:100];
    cfg.t_ftimwin    = ones(length(cfg.foi),1);
    cfg.toi          = 0:0.1:td.time{1}(end);
    ff               = ft_freqanalysis(cfg, td);  
    
    cfg             = [];
    cfg.method      = 'coh';
    coh             = ft_connectivityanalysis(cfg, ff); 

    cfg              = [];
    cfg.avgovertime  = 'yes';
    cfg.nanmean      = 'yes';
    coh              = ft_selectdata(cfg,coh);    
    
    out.Cohfrq        = coh.freq;
    out.Coh(cnt1,:)   = squeeze(coh.cohspctrm(1,2,:)).^2;
    
    cf                = (squeeze(out.SPLFP(1,:).^2)*Ntest-1)/(Ntest-1);
    out.CohAn(cnt1,:) = (cf*out.g(cnt1)*(Ntest-1)+out.g(cnt1))./(cf*(out.g(cnt1)*Ntest-1)+1);
    
    [val1,ind]          = max(squeeze(out.Coh(cnt1,:)),[],2); 
    out.CohPeak(cnt1)   = val1;
    out.CohPeakAn(cnt1) = out.CohAn(cnt1,ind);
end


end