% Model of Figure S6
% add fieldtrip toolbox to your path
clear all
addpath(fullfile(cd,'functions'));

%%
Test.N        = 1%500;
Test.cw       = [logspace(-3,0,15)];
Test.alpha    = [logspace(-2,0,10) 1.25 1.5 2];
Test.N        = 500;
Test.fac      = [0.242,0.39,0.556] ;

id = 1;
for cnt1= 1 : length(Test.cw)
    for cnt2 = 1 : length(Test.fac)
        cw{id}      = Test.cw(cnt1);
        alpha{id}   = Test.fac(cnt2); 
        g{id}       = 1;
        N{id}       = Test.N;
        steps{id}   = 8;
        ind{id}.sos = cnt2;
        ind{id}.cw  = cnt1;
        id          = id + 1;
    end
end

cellfun(@LIFmod,cw,alpha,g,N,steps,ind)
