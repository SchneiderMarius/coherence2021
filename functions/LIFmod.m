function [out] = LIFmod(par)

steps   = 10;   % different population size 
Ntest   = 80;   % set of neuron to evaluate PPC
Ncorr   = 0.5;  % factor scaling correlation
Nscale  = 13.3; % factor scaling LFP amp

steps = round(linspace(1,par.N,steps)/10)*10;
if length(steps)>1
    steps(1) = 1;
end

num     = 3;
delay	= 4;
time    = 1/par.fsample:1/par.fsample:par.time;
len     = length(time);

Noise       = cell(1,par.numTrial);
Signal      = cell(1,par.numTrial);
SPt         = cell(1,par.numTrial);
Frate       = zeros(1,par.numTrial);  
out.LFP     = cell(1,par.numTrial);
out.SP      = cell(length(steps),par.numTrial); 

for cnt0 = 1 : par.numTrial
    S = zeros(length(time),1);
    for cnt1 = 1 : size(par.coeff,1)
        S(1:2,cnt1) =  randn(2,1);
        for cnt2 = 3 : length(time)            
           S(cnt2,cnt1) = par.coeff(cnt1,:) * S(cnt2-2:cnt2-1,cnt1) + randn;     
        end
        pink = Pink_noise(len,num+par.N,par.a);
        
        Signal{cnt0} = par.fac*S;
        NoiseInit    = 1000*pink;
        
        out.LFP{cnt0}(:,1)          = Signal{cnt0}(:,1) + sqrt((1-par.g))*NoiseInit(:,1)+sqrt(par.g)*NoiseInit(:,2);
        LFP(:,2:par.N+1)            = sqrt(Ncorr)*NoiseInit(:,end)+sqrt(1-Ncorr)*NoiseInit(:,3:par.N+2);
        Noise{cnt0}                 = sum(LFP(:,2:par.N+1),2);        
        LFP(delay+1:end,2:par.N+1)  = sqrt(Ncorr)*NoiseInit(delay+1:end,end)+sqrt(1-Ncorr)*NoiseInit(delay+1:end,3:par.N+2) + par.cw*(Nscale*Signal{cnt0}(1:end-delay,1)/sqrt(par.N)); 
        Signal{cnt0}                = sum(Nscale*par.N*Signal{cnt0}/sqrt(par.N),2);          
    end

    input = (LFP(:,2:end)-mean(LFP(:,2:end),1))*0.038; %scales amplitude
        
    [~,spt,Frate_temp] = LIF(input,[],[],par.N);
    
    for cnt1 = 1 : length(steps)
        out.SP{cnt1,cnt0}       = zeros(1,size(LFP,1));
        out.SP{cnt1,cnt0}(1,:)  = full(sum(spt(:,1:steps(cnt1))',1)); 
    end
    Frate(cnt0) = mean(Frate_temp);
    SPt{cnt0} = full(spt(:,1:Ntest)');
end

clear pink S input Frate_temp spt
out.Frate  = mean(Frate);
out.params = par;
out.N      = steps;

for cnt1 = 1 : length(steps)

    [td] = fieldtrip(out.SP(cnt1,:),par);
    cfg          = [];
    cfg.lpfilter = 'yes';
    cfg.lpfreq   = 100;
    cfg.hpfilter = 'yes';
    cfg.hpfreq   = 1;
    [td]         = ft_preprocessing(cfg,td);
    
    for cntT = 1 : length(out.LFP)
        out.LFP{cntT}(:,2) = td.trial{cntT}(1,:);
    end  
	[out.f,~,~,out.Coh(cnt1,:),out.SOS,out.pow(:,:,cnt1),~,~] = AR_CohSOS(out.LFP,Signal,Noise,par,500,300);      
end
