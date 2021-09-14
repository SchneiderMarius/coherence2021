% SSM Model
clear all
addpath(fullfile(cd,'functions'));

load(fullfile(cd,'par','ARparameter'));

par             = parAR;
par.numTrial    = 500;
par.fsample     = 1000;
par.time        = 20;                                   % time in sec
par.a           = 2/3;                                  % slope 1/f

func    = {'sigmoid',[]};                               % determines transfer function
cw      = [logspace(-3,0,15) 2 4 6 8 10];               % projection strength
fac     = [logspace(-2,0,10) 2 4 6 8 10 14 50 100];     % determines SOS
result  = cell(1,length(func));

for cnt1 = 1 : length(func)
    
    id =1 ;
    parameter = cell(1,length(fac)*length(cw));
    for cnt2 = 1 : length(fac)        
        for cnt3 = 1 : length(cw)
            parameter{id}       = par;
            parameter{id}.cw    = cw(cnt3);
            parameter{id}.fac   = fac(cnt2);
            parameter{id}.func  = func{cnt1};
            id = id + 1;
        end
    end            
    out = cellfun(@SSM1,parameter);  
    
    id = 1;
    for cnt2 = 1 : length(fac)        
        for cnt3 = 1 : length(cw)
            result{cnt1}.cw(cnt2,cnt3) = cw(cnt3);
            result{cnt1}.sos(cnt2,cnt3) = out(id).SOS;
            result{cnt1}.CohSim(cnt2,cnt3) = out(id).Coh;  
            result{cnt1}.CohAna(cnt2,cnt3) = sqrt(1./(1+(1./(cw(cnt3).^2*(1+out(id).SOS)))));
            result{cnt1}.func = func{cnt1};
            id = id +1;
        end
    end
end

