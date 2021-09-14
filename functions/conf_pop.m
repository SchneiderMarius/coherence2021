function [W,init_state,input,alpha,beta,td] = conf_pop(par)


par.cell = [];
alpha = [];
beta  = [];
input = [];

par.synweight = zeros(length(par.celltype),length(par.celltype));

id = 1;
for cnt1 = 1 : length(par.regions)
    par.synweight(id:id+1,id:id+1)= [par.wEE(cnt1) par.wEI(cnt1); par.wIE(cnt1) par.wII(cnt1)];
    id = id + sum(contains(par.celltype,sprintf('%d',cnt1)));
end

for cnt1 = 1 : length(par.celltype)
    par.cell     = [par.cell; ones(par.cellnum(cnt1),1)*cnt1];
    alpha = [alpha; ones(par.cellnum(cnt1),1)* par.alpha(cnt1)];
    beta  = [beta; ones(par.cellnum(cnt1),1)* par.beta(cnt1)];
    input = [input; ones(par.cellnum(cnt1),1)* par.input(cnt1)];
    for cnt2 = 1 : length(par.celltype)
        if par.cellnum(cnt2)>0
            par.prob(cnt1,cnt2) = par.intraconn; 
            par.syn(cnt1,cnt2) = par.synweight(cnt1,cnt2)/(par.intraconn*par.cellnum(cnt2));    
        else
            par.prob(cnt1,cnt2) = 0; 
            par.syn(cnt1,cnt2) = 0;
        end
    end
end

for cnt1 = 1 : size(par.interconnlabel,1)
   ids = find(contains(par.celltype,par.interconnlabel{cnt1,1})); 
   idt = find(contains(par.celltype,par.interconnlabel{cnt1,2})); 

   
   if par.interconn(cnt1)>1
      par.syn(idt,ids)  = par.syn(ids,ids)*par.interconn(cnt1);
   else
      par.syn(idt,ids)  = par.syn(ids,ids);
   end

   par.prob(idt,ids) = par.interconn(cnt1); 
end


[W] = set_conn(par);

init_state(1,:) = round(rand(1,sum(par.cellnum)));
init_state(2,:) = abs(init_state(1,:)-1);

td = [];
td.fsample = par.fsample;

cnum = cumsum(par.cellnum);

for cnt = 1 : sum(par.cellnum)
    if cnt<=cnum(1)
       td.label{cnt} = sprintf('%s_%d',par.celltype{1},cnt); 
    else 
       id = find(cnum>=cnt,1);
       td.label{cnt} = sprintf('%s_%d',par.celltype{id},cnt-cnum(id-1)); 
    end
end

end