function [LFP,LFPsyn] = SP_LFP(SP,par) 

celltype=[0 cumsum(par.cellnum)];

cells   = cell(1,length(SP.trial));
for cnt1=1 : sum(par.cellnum>0)
    numT = celltype(cnt1)+1:celltype(cnt1+1); 
    numS = round(numel(celltype(cnt1)+1:celltype(cnt1+1)));
    out  = randperm(length(numT));
    cells{cnt1} = numT(out(1:numS));    
end


tau1 = 3;
tau2 = 2;
tsyn = 0 : 5*(tau1+tau2);
gsyn = ( exp(-tsyn/tau1)-exp(-tsyn/tau2));


LFPsyn      = cell(1,length(SP.trial));
activity    = cell(1,length(SP.trial));
LFP         = cell(1,length(SP.trial));

for cnt1 = 1 : length(SP.trial)
    
    LFPsyn{cnt1}           = zeros (sum(par.cellnum>0), size(SP.trial{cnt1},2)+length(tsyn));
    
    for cnt2 = 1 : sum(par.cellnum>0)
        if contains(par.celltype{cnt2},'E')
            fac = 1;
        elseif contains(par.celltype{cnt2},'I')
            fac = -1;
        end
       
        activity{cnt1}(cnt2,:) = full(sum(SP.trial{cnt1}(cells{cnt2},:),1));
        
        
        for cnt3 = 1 : size(activity{cnt1},2)       
            t1       = cnt3;
            tcut     = t1 + 5*(tau1+tau2);
            LFPsyn{cnt1}(cnt2, t1 : t1 + (tcut - t1)) = ...
                LFPsyn{cnt1}(cnt2, t1 : t1 + (tcut - t1)) + ...
                (fac.*gsyn (1 : tcut - t1 + 1)*activity{cnt1}(cnt2,cnt3));
        end
    end
    
    LFPsyn{cnt1} = LFPsyn{cnt1}(:,1:size(SP.trial{cnt1},2));
    
    for cnt2 = 1 : length(par.regions)
        idr = find(contains(par.celltype,sprintf('%d',cnt2)));  
        idp = find(contains(par.interconnlabel(:,2),sprintf('%d',par.regions(cnt2))));
        if any(par.cellnum(idr)>0)
            LFP{cnt1}(cnt2,:) = sum(LFPsyn{cnt1}(idr,:),1);
            for cnt3 = 1 : length(idp)
                ip = find(contains(par.celltype,par.interconnlabel{idp(cnt3),1}));
                if contains(par.interconnlabel{idp(cnt3),1},'E')
                    fac = 1;
                elseif contains(par.interconnlabel{idp(cnt3),1},'I')
                    fac = -1;
                end             
                LFP{cnt1}(cnt2,:) = LFP{cnt1}(cnt2,:) + [0 (LFPsyn{cnt1}(ip,1:end-1)*par.interconn(idp(cnt3))*fac)];
            end  
        end
    end
end

end