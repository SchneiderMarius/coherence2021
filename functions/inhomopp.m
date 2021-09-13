function [spikeMat, Frate,activity] = inhomopp(mod,modstr,par,maxfire,steps)

% mod: normalised modulating process 

% input: should have same length as mod and vary between 0 and 10

% max. firing rate

if nargin<2
    modstr = 1;
end
if nargin<3
   par.fsample = 1/1000;
end
if nargin<4
    maxfire = 2;
end

    
tstop = length(mod)/par.fsample;
t = 1/par.fsample : 1/par.fsample : length(mod)/par.fsample;
spikeMat = zeros(par.Ntest,length(t));

md = (modstr*mod' - mean(modstr*mod'));

for cnt1 = 1  : par.Ntest

    modt = md + maxfire;
    
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
Frate = sum(spikeMat,[1 2])/(par.Ntest*tstop);


for cnt1 = 1 : length(steps)
    activity(cnt1,:) = sum(spikeMat(1:steps(cnt1),:),1);
end
