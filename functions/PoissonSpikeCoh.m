function  [out] = PoissonSpikeCoh(par)

f = @(x) 2*std(x)*(2*x/std(x))./(1+2*abs(x/std(x)));   

num             = 3;
delay           = 4;
len             = par.fsample*par.time;
time            = 1/par.fsample:1/par.fsample:par.time;

par.steps       = round(linspace(1,par.N,par.stepnum)/10)*10;
par.steps(1)    = 1;

out.N = par.steps;

out.Signal  = cell(1,par.numTrial);
out.Noise   = cell(1,par.numTrial);
out.LFP     = cell(1,par.numTrial);

for cnt0 = 1 : par.numTrial
    S = zeros(length(time),1);
    for cnt1 = 1 : size(par.coeff,1)
        S(1:2,cnt1) =  randn(2,1);
        for cnt2 = 3 : length(time)            
           S(cnt2,cnt1) = par.coeff(cnt1,:) * S(cnt2-2:cnt2-1,cnt1) + randn;     
        end
        pink = Pink_noise(len,num,par.a);
        
        out.Signal{cnt0} = S;
        out.Noise{cnt0}  = par.fac*1000*pink;
        out.LFP{cnt0}(:,1) = out.Signal{cnt0}(:,1) + sqrt((1-par.g))*out.Noise{cnt0}(:,1)+sqrt(par.g)*out.Noise{cnt0}(:,2);
        out.LFP{cnt0}(:,2) = out.Noise{cnt0}(:,3);
    
        temp =  out.Noise{cnt0}(delay+1:end,3) + par.cw*(out.Signal{cnt0}(1:end-delay,1) + out.Noise{cnt0}(1:end-delay,1)); 
        out.LFP{cnt0}(delay+1:end,2) = f(temp);    
    end
end

clear S  pink time 

LF = cell(1,length(out.LFP));
for cnt1 = 1 : length(out.LFP)  
    [~, ~,LF{cnt1}]         = inhomopp(out.LFP{cnt1}(:,2),par);%)modstr,params,maxfire,steps);      
end

for cnt1 = 1 : length(par.steps)
    
    LFP = cell(1,length(LF));
    for cntT =1 : length(LF)
        LFP{cntT} = full(LF{cntT}(cnt1,:)); 
    end

    [td] = fieldtrip(LFP,par);
    cfg          = [];
    cfg.lpfilter = 'yes';
    cfg.lpfreq   = 100;
    cfg.hpfilter = 'yes';
    cfg.hpfreq   = 1;
    [td]         = ft_preprocessing(cfg,td);

    for cntT = 1 : length(out.LFP)
        out.LFP{cntT}(:,2) = td.trial{cntT}(1,:);
    end
    [out.f,~,~,out.Coh(cnt1,:),out.SOS,out.pow(:,:,cnt1),~,~] = AR_CohSOS(out.LFP,out.Signal,out.Noise,par,500,300);      
end

out = rmfield(out,'LFP');
out = rmfield(out,'Signal');
out = rmfield(out,'Noise');
end
