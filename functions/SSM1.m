function [out] =  SSM1(par)

if ~isfield(par,'func')
      f = @(x) x;       
else
    if strcmp(par.func,'sigmoid')
      f = @(x) min([max(x),abs(min(x))])*(x/std(x))./(1+abs(x/std(x)));
    else
      f = @(x) x;       
    end    
end

num         = 2;
delay       = 4;
len         = par.fsample*par.time;
time        = 1/par.fsample:1/par.fsample:par.time;
out.Signal  = cell(1,par.numTrial);
out.Noise   = cell(1,par.numTrial);
out.LFP     = cell(1,par.numTrial);

for cnt0 = 1 : par.numTrial
    S = zeros(length(time),1);
    for cnt1 = 1 : size(par.coeff,1)
        S(1:2,cnt1) = randn(2,1);
        for cnt2 = 3 : length(time)              
           S(cnt2,cnt1) = par.coeff(cnt1,:) * S(cnt2-2:cnt2-1,cnt1) + randn;         
        end
    end
    pink = Pink_noise(len,num,par.a);
    
    out.Signal{cnt0} = S;
    out.Noise{cnt0}  = par.fac*1000*pink;

    out.LFP{cnt0}(:,1) =  out.Signal{cnt0}(:,1) + out.Noise{cnt0}(:,1); 
    out.LFP{cnt0}(:,2) =  out.Noise{cnt0}(:,2); 
 
    temp = out.LFP{cnt0}(delay+1:end,2) + par.cw*out.LFP{cnt0}(1:end-delay,1);
    out.LFP{cnt0}(delay+1:end,2) = f(temp);
end

[~,~,~,out.Coh,out.SOS,~,~,~] = AR_CohSOS(out.LFP,out.Signal,out.Noise,par,500,300);      

out = rmfield(out,'LFP');
out = rmfield(out,'Signal');
out = rmfield(out,'Noise');
end   