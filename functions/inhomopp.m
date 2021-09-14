function [spikeMat, Frate,activity] = inhomopp(mod,par)
    
tstop       = length(mod)/par.fsample;
t           = 1/par.fsample : 1/par.fsample : length(mod)/par.fsample;
spikeMat    = zeros(par.N,length(t));

md = (par.modstr*mod' - mean(par.modstr*mod'));

for cnt1 = 1  : par.N

    modt = md + par.frate;
    
    m2 = max(modt);   
    u = rand(1,ceil(1.5*tstop*m2));
    y = cumsum(-(1/m2)*log(u));
    y = y(y<tstop); n=length(y);
    y1 = ceil(y*par.fsample);
    m = modt(y1);
    y = y(rand(1,n)<m/m2);
    spt = ceil(y*par.fsample);
    spikeMat(cnt1,spt) = 1;      
    spikeMat = sparse(spikeMat);
end
Frate = sum(spikeMat,[1 2])/(par.N*tstop);

activity = [];
for cnt1 = 1 : length(par.steps)
    activity(cnt1,:) = sum(spikeMat(1:par.steps(cnt1),:),1);
end
