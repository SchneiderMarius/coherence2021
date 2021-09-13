% SSM Model
clear all
addpath(fullfile(cd,'functions'));

%%
load(fullfile(cd,'par','ARparameter'));

numTrial        = 500;
num             = 2;
params.fsample  = 1000;
params.rate     = 1/params.fsample;
params.time     = 30;
params.delay    = 4;
params.a        = 2/3;
params.freq     = [17 19];
len             = 1/params.rate*params.time;
time            = params.rate:params.rate:params.time;

cw  = [logspace(-3,0,15) 2 4 6 8 10];
SNR = [logspace(-2,0,10) 2 4 6 8 10 14 50 100];

CohSim = zeros(length(SNR),length(cw),numTrial);
CohAna = zeros(length(SNR),length(cw),numTrial);

fac     = fliplr(SNR);
func    = {'sigmoid',[]};

result = cell(1,length(func));

for cnt1 = 1 : length(func)
    id =1 ;
    for cnt2 = 1 : length(SNR)        
        for cnt3 = 1 : length(cw)
            funcin{id}      = func{cnt1};
            parARin{id}     = parAR;
            paramsin{id}    = params;
            numTrialin{id}  = numTrial;
            facin{id}       = fac(cnt2);
            cwin{id}        = cw(cnt3);
            id = id + 1;
        end
    end            
    out = cellfun(@SSM,parARin,paramsin,numTrialin,facin,cwin,funcin)                 
    id = 1;
    for cnt2 = 1 : length(SNR)        
        for cnt3 = 1 : length(cw)
            result{cnt1}.cw(cnt2,cnt3) = cw(cnt3);
            result{cnt1}.snr(cnt2,cnt3) = out{id}.snr_temp;
            result{cnt1}.CohSim(cnt2,cnt3) = out{id}.Cohf;  
            result{cnt1}.cohAna(cnt2,cnt3) = sqrt(1./(1+(1./(cw(cnt3).^2*(1+out{id}.snr_temp)))));
            result{cnt1}.func = func{cnt1};
            id = id +1;
        end
    end
end
save(fullfile(cd,'results','Fig2BC'),'result');