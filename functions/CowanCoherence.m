function CowanCoherence(name,cntIn,cntOut,strength,apend)

data = dir(fullfile(cd,sprintf('%s*strngth%d_trl*',name,cntIn)));

for cnt0 = 1 : length(data)
    load(fullfile(data(cnt0).folder,data(cnt0).name));
    tdSP.trial{cnt0} = SPt;
    if cnt0 == 1
       param = par;
    end
end

param.interconn   = [strength];
input           = tdSP;


[LFP,~]     = SP_LFP(input,param); 
[~,tdLFP,~] = fieldtrip(LFP,param);

cfg              = [];
cfg.output       = 'pow';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = 1:1:100;                    
cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5; 
cfg.toi          = 0.1:0.05:3;         
freqP            = ft_freqanalysis(cfg, tdLFP);

cfg              = [];
cfg.output       = 'fourier';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = 1:1:100;                   
cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5; 
cfg.toi          = 0.1:0.05:3;                
freqF            = ft_freqanalysis(cfg, tdLFP);        

cfg = [];
cfg.avgovertime  = 'yes';
cfg.nanmean      = 'yes';
freqP            = ft_selectdata(cfg,freqP);

cfg              = [];
cfg.method       = 'coh';
coh              = ft_connectivityanalysis(cfg, freqF);        

cfg              = [];
cfg.avgovertime  = 'yes';
cfg.nanmean      = 'yes';
coh              = ft_selectdata(cfg,coh);  

result.pwrspctrm = freqP.powspctrm;
result.pwrfrq    = freqP.freq;
result.chspctrm  = coh.cohspctrm;
result.frq       = coh.freq;
result.strength  = strength;    

clear coh freqP

cfg              = [];
cfg.method       = 'granger';
granger          = ft_connectivityanalysis(cfg, freqF);        

cfg              = [];
cfg.avgovertime  = 'yes';
cfg.nanmean      = 'yes';
granger              = ft_selectdata(cfg,granger);  

result.grangerspctrm  = granger.grangerspctrm;
result.grangerfrq       = granger.freq;

save(fullfile(cd,sprintf('Coherence_%s_%s_strength%d',name,apend,cntOut)),'result');

end