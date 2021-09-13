function  [out] = PoissonSpikePPC(cw,alpha,maxfire,g,Ntest,steps)

if nargin<1
    cw = 0.05;
end
if nargin<2
   alpha = 0.9500; 
end
if nargin<3
    maxfire = 10000;
end
if nargin<4
    g = 0.95;    
end
if nargin<5
    Ntest = 50;    
end

if nargin<6
   steps = 1;
end
if nargin<7
   func = 'sigmoid';
end

load(fullfile(cd,'par','Fig2_ARfit'))

f = @(x) 2*std(x)*(x/std(x))./(1+abs(x/std(x)));


params.fsample  = 1000;
num             = 3;
params.rate     = 1/params.fsample;
params.time     = 10;
len             = 1/params.rate*params.time;
time            = params.rate:params.rate:params.time;
numTrial        = 1500;
delay = 4;

out.Noise       = cell(1,numTrial);
out.LFP         = cell(1,numTrial);

for cnt0 = 1 : numTrial
    S = zeros(length(time),1);
    for cnt1 = 1 : size(parAR{1}.coeff,1)
        S(1:2,cnt1) =  randn(2,1);
        for cnt2 = 3 : length(time)            
           S(cnt2,cnt1) = parAR{1}.coeff(cnt1,:) * S(cnt2-2:cnt2-1,cnt1) + randn;     
        end
        pink = Pink_noise(len,num,parAR{1}.a,params);
        
        out.Signal{cnt0} = S;
        out.Noise{cnt0}  = alpha*1000*pink;
        out.LFP{cnt0}(:,1) = out.Signal{cnt0}(:,1) + sqrt((1-g))*out.Noise{cnt0}(:,1)+sqrt(g)*out.Noise{cnt0}(:,2);
        out.LFP{cnt0}(:,2) = out.Noise{cnt0}(:,3);
    
        temp =  out.Noise{cnt0}(delay+1:end,3) + cw*(out.Signal{cnt0}(1:end-delay,1) + out.Noise{cnt0}(1:end-delay,1)); 
        out.LFP{cnt0}(delay+1:end,2) = f(temp);    
    end
end

[~,~,~,~,out.SNR,~,~,~] = AR_CohSOS(out.LFP,out.Signal,out.Noise,[],params,500,300);      

modstr      = 0.1;
steps       = round(linspace(1,Ntest,steps+1)/10)*10;
steps(1)    = 1;

params.Ntest = Ntest;
SP = cell(1,length(out.LFP)); LF = cell(1,length(out.LFP));
for cnt1 = 1 : length(out.LFP)  
    [spikeMat, Frate(cnt1),activity]  = inhomopp(out.LFP{cnt1}(:,2),modstr,params,maxfire,steps);      
    SP{cnt1}                          = sum(full(spikeMat(1:Ntest,:)),1);
    LF{cnt1}                          = full(activity);
end

%%
idlfp{1} = 1;
idlfp{2} = 2;
idN{1}   = 3;

SPLFP = [];
for cntT = 1 : length(out.LFP)
    SPLFP{cntT} = [out.LFP{cntT}(:,1)';out.LFP{cntT}(:,2)';SP{cntT}];
end

[~,tdSP,~] = fieldtrip(SPLFP,params);

clear td SPLFP S SP LF input1 pink Projection spikeMat time activity
out = rmfield(out,'LFP');
out = rmfield(out,'Signal');
out = rmfield(out,'Noise');

out.PPC = [];
out.Coh = [];

idc = 1;
for cntL = 1 : length(idlfp)
    for cntS = 1 : length(idN)
        PPC_temp = [];
        Coh_temp = [];
        for cnt2 = 1 : length(idN{cntS})
            cfg              = [];
            cfg.method       = 'mtmfft';
            cfg.foi          = [1:100];
            cfg.taper        = 'hanning'; 
            cfg.timwin       = [-0.35/2 0.35/2];
            cfg.spikechannel = tdSP.label(idN{cntS}(cnt2));
            cfg.channel      = tdSP.label(idlfp{cntL});
            stsFFT           = ft_spiketriggeredspectrum(cfg, tdSP) 
            
            cfg               = [];
            cfg.method        = 'ppc0';
            cfg.spikechannel  = tdSP.label(idN{cntS}(cnt2));
            cfg.channel       = tdSP.label(idlfp{cntL});
            cfg.avgoverchan   = 'unweighted';
            statSts1          = ft_spiketriggeredspectrum_stat(cfg,stsFFT);   

            PPC_temp(cnt2,:)  = statSts1.ppc0;    
            out.PPCfrq = statSts1.freq;

            clear stsFFT statSts1
            
            cfg              = [];
            cfg.channel      = tdSP.label([idlfp{cntL} idN{cntS}(cnt2)]);
            cfg.output       = 'fourier';
            cfg.method       = 'mtmfft';
            cfg.taper        = 'hanning';
            cfg.foi          = [1:1:100]; 
            cfg.t_ftimwin    = 0.35;
            cfg.toi          = 0:0.05:tdSP.time{1}(end);
            fourier1         = ft_freqanalysis(cfg, tdSP);  
            
            cfg             = [];
            cfg.method      = 'coh';
            conn            = ft_connectivityanalysis(cfg, fourier1); 
            
            Coh_temp(cnt2,:)  = conn.cohspctrm(1,2,:);                
            out.Cohfrq        = conn.freq;
            clear fourier1 conn  
        end
        out.PPC(idc,:) = mean(PPC_temp,1);
        out.Coh(idc,:) = mean(Coh_temp,1);
        idc = idc + 1;
    end
end
out.Frate  = mean(Frate);
end
