% Model of Figure S8,S2,S5
% add fieldtrip toolbox to your path
clear all
addpath(fullfile(cd,'functions'));


%%
data{1} = dir(fullfile(cd,'par','*25Hz.mat'));
data{2} = dir(fullfile(cd,'par','*25Hz.mat'));
data{3} = dir(fullfile(cd,'par','*65Hz.mat'));

numTrial = 800;
strength = [0 logspace(-2,0,10) 0.1];
strength = sort(strength);
numSim   = length(strength)*numTrial;

for cnt1 = 1 : length(data)
    load(fullfile(data{cnt1}(1).folder,data{cnt1}(1).name));

    if cnt1~=2
        par.interconnlabel = {'E1','E2'};        
        name = sprintf('65_to_%d',par.name);
        result.snd = 1;
        result.rcv = 2;         
    elseif cnt1==2
        par.interconnlabel = {'E2','E1'};
        name = sprintf('%d_to_65',par.name);
        result.snd = 2;
        result.rcv = 1;         
    end
    
    
    param       = cell(1,numSim);
    W           = cell(1,numSim);
    init_state  = cell(1,numSim);
    input       = cell(1,numSim);
    alpha       = cell(1,numSim);
    beta        = cell(1,numSim);
    strngth     = cell(1,numSim);
    flder       = cell(1,numSim);
    response_fn = cell(1,numSim);
    id = 1;
    for cnt2 = 1 : length(strength)
        for cnt3 = 1: numTrial
            param{id} = par;
            param{id}.interconn = [strength(cnt2)];
            [W{id},init_state{id},input{id},alpha{id},beta{id},~] = conf_pop(param{id});
            strngth{id} = strength(cnt2);
            flder{id} = fullfile(cd,sprintf('%s_strngth%d_trl%d',name,cnt2,cnt3));
            response_fn{id}     = 'sigmoid';
            id = id + 1;
        end
    end
	cellfun(@WilsonCowan,W, input, alpha, beta,init_state, param,response_fn,flder);            
    
    
    apend   = cell(1,2*length(strength)); 
    strg    = cell(1,2*length(strength)); 
    cntOut  = cell(1,2*length(strength)); 
    cntIn   = cell(1,2*length(strength)); 
    nam     = cell(1,2*length(strength)); 
    
    for cnt2 = 1 : length(strength)
        nam{cnt2}  = name;
        cntIn{cnt2}  = cnt2;
        cntOut{cnt2}  = cnt2;
        strg{cnt2} = strength(cnt2);
        apend{cnt2} = 'entrain';
    end

    ct = length(strength);

    for cnt2 = 1 : length(strength)
        nam{cnt2+ct}  = name;
        cntIn{cnt2+ct}  = 1;
        cntOut{cnt2+ct}  = cnt2;
        strg{cnt2+ct} = strength(cnt2);
        apend{cnt2+ct} = 'unentrain';
    end
    
    cellfun(@CowanCoherence,nam,cntIn,cntOut,strg,apend)  
end
