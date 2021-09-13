function [f1,td,dat] = fieldtrip(data,params,pl)

if ~isfield(params,'granger')
    params.granger = 0;
end
if ~isfield(params,'coh')
    params.coh = 0;
end
if ~isfield(params,'avg')
    params.avg = 0;
end
if ~isfield(params,'pow')
    params.pow = 0;
end
if ~isfield(params, 'psd')
   params.psd = 0; 
end
if ~isfield(params, 'transfer')
   params.transfer = 0; 
end

dat=[];
f1 = [];

if ~isfield(params,'ft') || params.ft~=1 

    td = [];
    td.fsample = params.fsample;

    if ~iscell(data)
        time = 1/td.fsample:1/td.fsample:(size(data,1)-1)/td.fsample;
        td.time{1} = time;     
        td.trial{1} = data';
        td.sampleinfo(1,:) = [1 length(td.time{1})];
        for cnt = 1 : size(data,2)
            td.label{cnt} = sprintf('CH%d',cnt);
        end   
    else
        if size(data{1},1)>size(data{1},2)
            time = 1/td.fsample:1/td.fsample:size(data{1},1)*1/td.fsample;
        else
            time = 1/td.fsample:1/td.fsample:size(data{1},2)*1/td.fsample; 
        end
        for cnt = 1 : length(data)
            td.time{cnt} = time;
            if size(data{cnt},2)==length(time)
                td.trial{cnt} = data{cnt};
            else
                td.trial{cnt} = data{cnt}';            
            end
            if cnt>1
                td.sampleinfo(cnt,:) = [td.sampleinfo(cnt-1,2)+1 length(td.time{1})*cnt];
            else
                td.sampleinfo(cnt,:) = [1 length(td.time{1})];            
            end
        end
        for cnt = 1 : size(td.trial{1},1)
            td.label{cnt,1} = sprintf('CH%d',cnt);
        end       
    end

elseif params.ft==1  
    td = data;
end


if nargin<3
    f1=[];
end

if nargin==3 || params.pow==1

    frq             = 1 : 1 : 100;
    cfg              = [];
    cfg.output       = 'pow';  
    cfg.method       = 'mtmfft';
    cfg.foi          = frq; 
    cfg.taper        = 'hanning'; 
    cfg.t_ftimwin    = ones(length(cfg.foi),1)*0.35;
    cfg.toi          = 0.1:0.05:td.time{1}(end); 
    dat.freq         = ft_freqanalysis(cfg, td); 

    
    if params.avg == 1 && isfield(dat.freq,'time')
        cfg             = [];
        cfg.avgovertime = 'yes';
        cfg.nanmean     = 'yes';
        dat.freq        = ft_selectdata(cfg, dat.freq)
    end
end

if params.psd==1
    freqP=[]; 
    frq             = 1 : 1 : 100;
    cfg             = [];
    cfg.output      = 'fourier';
    cfg.method      = 'mtmconvol';
    cfg.taper       = 'hanning';
    cfg.foi         = frq; 
    cfg.tapsmofrq   = 2;
    numfoi          = length(cfg.foi);
    cfg.t_ftimwin   = 0.5*ones(numfoi);
    cfg.toi         = 0:1:0.05:time(end);
    dat.freqP       = ft_freqanalysis(cfg, td);
    
    dat.freqP.fourierspctrm = dat.freqP.fourierspctrm.*conj(dat.freqP.fourierspctrm);

    cfg = []
    cfg.avgovertime = 'yes'
    cfg.nanmean = 'yes';
    cfg.avgoverrpt = 'yes'
    dat.freqP = ft_selectdata(cfg, dat.freqP)
end


if nargin==3 && any(contains(pl,'freq'))
   f1 = figure
   plot(freqP.freq,freqP.powspctrm)
   xlabel('frequency [Hz]')
   ylabel('power')
elseif nargin==3 && (any(contains(pl,'coh')) || any(contains(pl,'granger'))) || params.coh ==1

    frq              = 1 : 1 : 100;
    cfg              = [];
    cfg.output       = 'fourier';
    cfg.method       = 'mtmconvol';
    cfg.taper        = 'hanning';
    cfg.foi          = frq;
    cfg.t_ftimwin    = ones(length(cfg.foi),1)*0.35;
    cfg.toi          = 0.1:0.05:time(end);    
    freq             = ft_freqanalysis(cfg, td); 
    cfg             = [];
    cfg.method      = 'coh';
    dat.coh         = ft_connectivityanalysis(cfg, freq);
    
    if params.avg == 1 && isfield(dat.coh,'time')
        cfg = [];
        cfg.avgovertime = 'yes';
        cfg.nanmean = 'yes';
        dat.coh = ft_selectdata(cfg, dat.coh);
    end
end

if nargin==3 && any(strcmp(pl,'granger')) || params.granger==1
    cfg           = [];
    cfg.method    = 'mtmfft';
    cfg.taper     = 'dpss';
    cfg.output    = 'fourier';
    cfg.foilim    = [0 150];
    cfg.tapsmofrq = 2;
    freq          = ft_freqanalysis(cfg, td);    
    cfg         = [];
    cfg.method  = 'granger';
    dat.granger     = ft_connectivityanalysis(cfg, freq);

  
    if params.transfer == 1
        cfg                     = [];
        cfg.method              = 'granger';
        cfg.marius              = 3;
        dat.transfer            = ft_transfer(cfg, freq);
        

        siz(4:5) = [size(dat.transfer.crsspctrm,3) size(dat.transfer.crsspctrm,4)];
        for ii = 1 : size(dat.transfer.crsspctrm,1)
            for jj = 1 : size(dat.transfer.crsspctrm,2)
                Z  = dat.transfer.noisecov; 
                zc = reshape(Z(jj,jj,:) - Z(ii,jj,:).^2./Z(ii,ii,:),[1 1 1 siz(5)]);
                zc = repmat(zc,[1 1 siz(4) 1]);
                S  = dat.transfer.crsspctrm;
                H  = dat.transfer.transfer;
                Syy  = reshape(abs(S(ii,ii,:,:)),[1 1 siz(4:end)]);
                Sxx  = reshape(abs(S(jj,jj,:,:)),[1 1 siz(4:end)]);
                Syx  = reshape(abs(zc.*abs(H(ii,jj,:,:)).^2),[1 1 siz(4:end)]);    
                Sxy  = reshape(abs(zc.*abs(H(ii,jj,:,:)).^2),[1 1 siz(4:end)]);    
                transfer_temp = (Syy./(Syy-Syx)).*(Syx./Sxx);
                dat.ptransfer(jj,ii,:) = squeeze(nanmean(transfer_temp,4));   
                dat.S(jj,ii,:) = Syx;
                dat.S(ii,ii,:) = Syy;
                dat.S(jj,jj,:) = Sxx;
            end
        end
    end
elseif nargin==3 && any(strcmp(pl,'model_granger')) && params.plot==1
    cfg         = [];
    cfg.order   = 5;
    cfg.toolbox = 'bsmart';
    mdata       = ft_mvaranalysis(cfg, td);
    cfg         = [];
    cfg.method  = 'mvar';
    cfg.foi     = [1:100]    
    mfreq       = ft_freqanalysis(cfg, mdata);        
    cfg         = [];
    cfg.method  = 'granger';
    granger     = ft_connectivityanalysis(cfg, mfreq);
    f1 = figure
    subplot(3,1,1);  plot(freqP.freq,freqP.powspctrm); xlabel('freqenzy [Hz]'); ylabel('power')
    subplot(3,1,2); plot(conn.freq,squeeze(conn.cohspctrm(1,2,:))); xlabel('frequency [Hz]'); ylabel('coherence')
    subplot(3,1,3); plot(granger.freq,squeeze(granger.grangerspctrm(1,2,:))); xlabel('frequency [Hz]'); ylabel('granger')
    subplot(3,1,3); plot(granger.freq,squeeze(granger.grangerspctrm(2,1,:))); xlabel('frequency [Hz]'); ylabel('granger')

elseif nargin==3 && contains(pl,'coh') && params.plot==1
    f1 = figure
    subplot(2,1,1);  plot(freqP.freq,freqP.powspctrm); xlabel('freqenzy [Hz]'); ylabel('power')
    subplot(2,1,2); plot(conn.freq,squeeze(conn.cohspctrm(1,2,:))); xlabel('frequency [Hz]'); ylabel('coherence')
end
    
end
