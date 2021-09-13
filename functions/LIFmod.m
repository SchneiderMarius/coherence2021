function LIFmod(cw,alpha,g,N,steps,ind,numTrial)

if nargin<1
   cw = 0.05;
end
if nargin<2
   alpha = 0.9500; 
end
if nargin<3 || isempty(g)
    g = 0.95;    
end
if nargin<4 || isempty(N)
    N = 100;
end
if nargin<5 || isempty(N)
    steps = 1;
end
if nargin<7
   numTrial = 800;
end

load(fullfile(cd,'par','ARparameter'));

Ntest   = 80;
Ncorr   = 0.5;
Nscale = 13.3;

steps = round(linspace(1,N,steps)/10)*10;
if length(steps)>1
    steps(1) = 1;
end

params.fsample  = 1000;
params.time     = 4;
params.rate     = 1/params.fsample;
delay           = 4;

num     = 3;
len     = round(1/params.rate*params.time);
time    = params.rate:params.rate:params.time;
Signal  = cell(1,numTrial);
Noise   = cell(1,numTrial);
out.LFP = cell(1,numTrial);
out.SP  = cell(length(steps),numTrial);


LFPtest     = cell(1,numTrial);
allN        = cell(1,numTrial);
singleN     = cell(1,numTrial);
SPt         = cell(1,numTrial);
Frate       = zeros(1,numTrial);  
    
for cnt0 = 1 : numTrial
    S = zeros(length(time),1);
    for cnt1 = 1 : size(parAR.coeff,1)
        S(1:2,cnt1) =  randn(2,1);
        for cnt2 = 3 : length(time)            
           S(cnt2,cnt1) = parAR.coeff(cnt1,:) * S(cnt2-2:cnt2-1,cnt1) + randn;     
        end
        pink = Pink_noise(len,num+N,parAR.a);
        
        Signal{cnt0} = alpha*S;
        NoiseInit  = 1000*pink;
        out.LFP{cnt0}(:,1) = Signal{cnt0}(:,1) + sqrt((1-g))*NoiseInit(:,1)+sqrt(g)*NoiseInit(:,2);
        
        LFP(:,2:N+1) = sqrt(Ncorr)*NoiseInit(:,end)+sqrt(1-Ncorr)*NoiseInit(:,3:N+2);
        
        LFPtest{cnt0} = LFP(:,2:3);
        Noise{cnt0}  = sum(LFP(:,2:N+1),2);        
        
        LFP(delay+1:end,2:N+1) =sqrt(Ncorr)*NoiseInit(delay+1:end,end)+sqrt(1-Ncorr)*NoiseInit(delay+1:end,3:N+2) + cw*(Nscale*Signal{cnt0}(1:end-delay,1)/sqrt(N)); 

        Signal{cnt0} = sum(Nscale*N*Signal{cnt0}/sqrt(N),2);

        allN{cnt0}    = sum(LFP(:,2:N+1),2);
        singleN{cnt0} = LFP(:,2);          
    end


     input = (LFP(:,2:end)-mean(LFP(:,2:end),1))*0.038;
        
    [~,spt,Frate_temp] = LIF(input,[],[],N);
    
    for cnt1 = 1 : length(steps)
        out.SP{cnt1,cnt0}       = zeros(1,size(LFP,1));
        out.SP{cnt1,cnt0}(1,:)  = full(sum(spt(:,1:steps(cnt1))',1)); 
    end
    Frate(cnt0) = mean(Frate_temp);
    SPt{cnt0} = full(spt(:,1:Ntest)');
end

[out.foutN,~,~,out.powAllN,out.powSingleN] = AR_SOS(allN,singleN,[],params,500,300);
[out.foutSNR,~,out.SNR,~,~]                = AR_SOS(Signal,Noise,[],params,500,300);

params.freq = 19;
[out.foutRec,out.CohRec,~,out.powRec] = AR_Coh(LFPtest,[],params,500,300);

clear pink S input Frate_temp spt V Signal Noise LFPtest
out.Frate  = mean(Frate);
out.params = params;
out.Frate  = mean(Frate);
out.N      = steps;

for cnt1 = 1 : length(steps)

    [~,td,~] = fieldtrip(out.SP(cnt1,:),params);
    cfg          = [];
    cfg.lpfilter = 'yes';
    cfg.lpfreq   = 100;
    cfg.hpfilter = 'yes';
    cfg.hpfreq   = 1;
    [td]         = ft_preprocessing(cfg,td);
    
    for cntT = 1 : length(out.LFP)
        out.LFP{cntT}(:,2) = td.trial{cntT}(1,:);
    end
    params.freq = 19;
    [out.foutCoh,out.Coh(:,cnt1),out.Cohf(cnt1),out.pow(:,:,cnt1)] = AR_Coh(out.LFP,[],params,500,300);      
end

tdSP = cell(1,numTrial);
for cnt1 = 1 : numTrial
    tdSP{cnt1} = [out.LFP{cnt1}';SPt{cnt1}];
end

[~,tdSP,~] = fieldtrip(tdSP,params);

idlfp = [1 2];
idsp = 3:Ntest+2;

out.SPLFPcoh = [];
for cnt1= 1 : length(idlfp)

    chn = cell(length(idsp),2);
    for cnt2 = 1 : length(idsp)
        chn{cnt2,1} = tdSP.label{idlfp(cnt1)};
        chn{cnt2,2} = tdSP.label{idsp(cnt2)};  
    end   
    coh_temp = [];
    for cnt2 = 1 : length(idsp)
        cfg              = [];
        cfg.output       = 'fourier';
        cfg.method       = 'mtmconvol';
        cfg.taper        = 'hanning';
        cfg.foi          =  [1:100];
        cfg.t_ftimwin    = ones(length(cfg.foi),1); 
        cfg.toi          = 0:0.1:tdSP.time{1}(end); 
        cfg.channel      = chn(cnt2,:);   
        ff               = ft_freqanalysis(cfg, tdSP); 
        
        cfg             = [];
        cfg.channelcmb  = chn(cnt2,:);
        cfg.method      = 'coh';
        coh             = ft_connectivityanalysis(cfg, ff); 

        cfg              = [];
        cfg.avgovertime  = 'yes';
        cfg.nanmean      = 'yes';
        coh              = ft_selectdata(cfg,coh);
        
        coh_temp(cnt2,:) = squeeze(coh.cohspctrm(1,2,:));
    end
    out.SPLFPcoh(cnt1,:) = nanmean(coh_temp,1);
end
out.SPLFPfreq = coh.freq;
clear coh ff

tdSP = [];
for cnt1 = 1 : numTrial
    tdSP{cnt1} = [out.LFP{cnt1}';sum(SPt{cnt1},1)];
end

out = rmfield(out,'LFP');

[~,tdSP,~] = fieldtrip(tdSP,params);

idlfp = [1 2];
idsp = 3;

for cnt1= 1 : length(idlfp)
    PPC_temp = [];
    for cnt2 = 1 : length(idsp)
        cfg              = [];
        cfg.method       = 'mtmfft';
        cfg.foi          = [1:100]; 
        cfg.taper        = 'hanning'; 
        cfg.timwin       = [-0.35/2 0.35/2];
        cfg.spikechannel = tdSP.label{idsp(cnt2)};
        cfg.channel      = tdSP.label{idlfp(cnt1)};
        stsFFT           = ft_spiketriggeredspectrum(cfg, tdSP);

        cfg               = [];
        cfg.method        = 'ppc1';
        cfg.spikechannel = tdSP.label{idsp(cnt2)};
        cfg.channel      = tdSP.label{idlfp(cnt1)};
        cfg.avgoverchan   = 'unweighted';
        statSts1          = ft_spiketriggeredspectrum_stat(cfg,stsFFT);   
        PPC_temp(cnt2,:)  = statSts1.ppc1;  
    end
    out.PPC(cnt1,:) = mean(PPC_temp,1);
end
out.PPCfreq(1,:) = statSts1.freq;

par.cw       = cw;
par.g        = g;
par.alpha    = alpha;
par.N    	 = N;
par.steps    = steps;
save(fullfile(cd,sprintf('LIF_N%d_cw%d_sos%d',N,ind.cw,ind.sos)),'par','out');
end