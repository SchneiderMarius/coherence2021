function  [out] = SpikeModPartTransferCoh(cw,alpha,maxfire,g,Ntest,steps,folder)

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
   steps = 10;
end
if nargin<7
   func = 'sigmoid';
end

load(fullfile(cd,'par','ARparameter'));
f = @(x) 2*std(x)*(2*x/std(x))./(1+2*abs(x/std(x)));   

params.fsample  = 1000;
num             = 3;
params.rate     = 1/params.fsample;
params.time     = 10;
len             = 1/params.rate*params.time;
time            = params.rate:params.rate:params.time;
numTrial        = 500;
delay           = 4;

out.Signal  = cell(1,numTrial);
out.Noise   = cell(1,numTrial);
out.LFP     = cell(1,numTrial);

for cnt0 = 1 : numTrial
    S = zeros(length(time),1);
    for cnt1 = 1 : size(parAR.coeff,1)
        S(1:2,cnt1) =  randn(2,1);
        for cnt2 = 3 : length(time)            
           S(cnt2,cnt1) = parAR.coeff(cnt1,:) * S(cnt2-2:cnt2-1,cnt1) + randn;     
        end
        pink = Pink_noise(len,num,parAR.a);
        
        out.Signal{cnt0} = S;
        out.Noise{cnt0}  = alpha*1000*pink;
        out.LFP{cnt0}(:,1) = out.Signal{cnt0}(:,1) + sqrt((1-g))*out.Noise{cnt0}(:,1)+sqrt(g)*out.Noise{cnt0}(:,2);
        out.LFP{cnt0}(:,2) = out.Noise{cnt0}(:,3);
    
        temp =  out.Noise{cnt0}(delay+1:end,3) + cw*(out.Signal{cnt0}(1:end-delay,1) + out.Noise{cnt0}(1:end-delay,1)); 
        out.LFP{cnt0}(delay+1:end,2) = f(temp);    
    end
end

%% 
clear S  input1 pink Projection time activity

modstr      = 0.1;
steps       = round(linspace(1,Ntest,steps)/10)*10;
steps(1)    = 1;

%%
params.Ntest = Ntest;
LF = cell(1,length(out.LFP));
for cnt1 = 1 : length(out.LFP)  
    input1                            = zeros(length(out.LFP{cnt1}),1);
    [~, Frate(cnt1),LF{cnt1}]         = inhomopp(out.LFP{cnt1}(:,2),modstr,params,maxfire,steps);      
end

%%
for cnt1 = 1 : length(steps)
    for cntT =1 : length(LF)
        LFP{cntT} = full(LF{cntT}(cnt1,:)); 
    end

    [~,td,~] = fieldtrip_it(LFP,params);
    cfg          = [];
    cfg.lpfilter = 'yes';
    cfg.lpfreq   = 100;
    cfg.hpfilter = 'yes';
    cfg.hpfreq   = 1;
    [td]         = ft_preprocessing(cfg,td);


    for cntT = 1 : length(out.LFP)
        out.LFP{cntT}(:,2) = td.trial{cntT}(1,:);
    end
    [out.f,~,~,out.Cohf(cnt1,:),out.SNR,out.pow(:,:,cnt1),~,~] = AR_CohSOS(out.LFP,out.Signal,out.Noise,[],params,500,300);      
end

clear td  S  LF input1 pink Projection time activity

out.Frate  = mean(Frate);
out.N      = steps;
out = rmfield(out,'LFP');
out = rmfield(out,'Signal');
out = rmfield(out,'Noise');
end
