function [out] =  SSM(parAR,params,numTrial,fac,cw,func)

if strcmp(func,'sigmoid')
   f = @(x) min([max(x),abs(min(x))])*(x/std(x))./(1+abs(x/std(x)));
else
   f = @(x) x;       
end

num         = 2;
params.time = 20;
len         = 1/params.rate*params.time;
time        = params.rate:params.rate:params.time;

out.Signal  = cell(1,numTrial);
out.Noise   = cell(1,numTrial);
out.LFP     = cell(1,numTrial);

for cnt0 = 1 : numTrial
    S = zeros(length(time),1);
    for cnt1 = 1 : size(parAR.coeff,1)
        S(1:2,cnt1) = randn(2,1);
        for cnt2 = 3 : length(time)              
           S(cnt2,cnt1) = parAR.coeff(cnt1,:) * S(cnt2-2:cnt2-1,cnt1) + randn;         
        end
    end
    pink = Pink_noise(len,num,parAR.a);
    
    delID = round(params.delay/1000/params.rate);

    out.Signal{cnt0} = S;
    out.Noise{cnt0}  = fac*1000*pink;

    out.LFP{cnt0}(:,1) =  out.Signal{cnt0}(:,1) + out.Noise{cnt0}(:,1); 
    out.LFP{cnt0}(:,2) =  out.Noise{cnt0}(:,2); 
 
    temp = out.LFP{cnt0}(delID+1:end,2) + cw*out.LFP{cnt0}(1:end-delID,1);
    out.LFP{cnt0}(delID+1:end,2) = f(temp);
end
[~,~,~,out.Cohf,out.snr_temp,~,~,~] = AR_CohSOS(out.LFP,out.Signal,out.Noise,[],params,500,300);      

out = rmfield(out,'LFP');
out = rmfield(out,'Signal');
out = rmfield(out,'Noise');
end   