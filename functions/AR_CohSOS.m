function [fout,Coh,SNR,Cohf,SNRf,pow,S,Ns] = AR_CohSOS(LFP,Signal,Noise,foi,params,window,ovl)

if nargin<6
    window = 1000;
    ovl = window/2;
elseif nargin<7
    ovl = window/2;    
end

coh_temp = [];  

if iscell(LFP)
    for cnt1 = 1 : length(LFP)
      x{cnt1}    = LFP{cnt1}(1:end,1);
      y{cnt1}    = LFP{cnt1}(1:end,2);
      Signal{cnt1} = Signal{cnt1}(1:end,1);
      Noise{cnt1}  = Noise{cnt1}(1:end,1); 

      z        = buffer(1:numel(x{cnt1}),window,ovl);
      xt{cnt1} = x{cnt1}(z(:,all(z)));   

      z        = buffer(1:numel(y{cnt1}),window,ovl);
      yt{cnt1} = y{cnt1}(z(:,all(z)));     

      z        = buffer(1:numel(Signal{cnt1}),window,ovl);
      St{cnt1} = Signal{cnt1}(z(:,all(z)));  

      z        = buffer(1:numel(Noise{cnt1}),window,ovl);
      Nt{cnt1} = Noise{cnt1}(z(:,all(z)));    

      N = size(Nt{cnt1},2);

      Xf{cnt1} = fft(xt{cnt1},window,1);
      Yf{cnt1} = fft(yt{cnt1},window,1);   

      Sf{cnt1} = fft(St{cnt1},window,1);
      Nf{cnt1} = fft(Nt{cnt1},window,1);    

      Xf{cnt1} = Xf{cnt1}/window;
      Yf{cnt1} = Yf{cnt1}/window;
      Sf{cnt1} = Sf{cnt1}/window;
      Nf{cnt1} = Nf{cnt1}/window;
      
      Xf{cnt1} = Xf{cnt1}(1:window/2+1,:);
      Yf{cnt1} = Yf{cnt1}(1:window/2+1,:);
      Sf{cnt1} = Sf{cnt1}(1:window/2+1,:);
      Nf{cnt1} = Nf{cnt1}(1:window/2+1,:);
      
      Xf{cnt1}(2:end-1,:) = 2*Xf{cnt1}(2:end-1,:);
      Yf{cnt1}(2:end-1,:) = 2*Yf{cnt1}(2:end-1,:);
      Sf{cnt1}(2:end-1,:) = 2*Sf{cnt1}(2:end-1,:);
      Nf{cnt1}(2:end-1,:) = 2*Nf{cnt1}(2:end-1,:);
          

      xy(cnt1,:) = sum(Xf{cnt1} .* conj(Yf{cnt1}), 2)/N;
      xx(cnt1,:) = sum(Xf{cnt1} .* conj(Xf{cnt1}), 2)/N;
      yy(cnt1,:) = sum(Yf{cnt1} .* conj(Yf{cnt1}), 2)/N;
      ss(cnt1,:) = sum(Sf{cnt1} .* conj(Sf{cnt1}), 2)/N;
      nn(cnt1,:) = sum(Nf{cnt1} .* conj(Nf{cnt1}), 2)/N;
        
      Pow_temp1(cnt1,:) = sum(abs(Xf{cnt1}),2)/N;
      Pow_temp2(cnt1,:) = sum(abs(Yf{cnt1}),2)/N;
      S_temp(cnt1,:)    = sum(abs(Sf{cnt1}),2)/N;
      N_temp(cnt1,:)    = sum(abs(Nf{cnt1}),2)/N;  

    end
elseif ~iscell(LFP)
      cnt1 = 1;  
      x{cnt1}    = LFP(101:end,1);
      y{cnt1}    = LFP(101:end,2);
      Signal2{cnt1} = Signal(101:end,1);
      Noise2{cnt1}  = Noise(101:end,1);

      z        = buffer(1:numel(x{cnt1}),window,ovl);
      xt{cnt1} = x{cnt1}(z(:,all(z)));   

      z        = buffer(1:numel(y{cnt1}),window,ovl);
      yt{cnt1} = y{cnt1}(z(:,all(z)));     

      z        = buffer(1:numel(Signal2{cnt1}),window,ovl);
      St{cnt1} = Signal2{cnt1}(z(:,all(z)));  

      z        = buffer(1:numel(Noise2{cnt1}),window,ovl);
      Nt{cnt1} = Noise2{cnt1}(z(:,all(z)));    

      N = size(Nt{cnt1},2);

      Xf{cnt1} = fft(xt{cnt1},window,1);
      Yf{cnt1} = fft(yt{cnt1},window,1);   

      Sf{cnt1} = fft(St{cnt1},window,1);
      Nf{cnt1} = fft(Nt{cnt1},window,1);    

      Xf{cnt1} = Xf{cnt1}/window;
      Yf{cnt1} = Yf{cnt1}/window;
      Sf{cnt1} = Sf{cnt1}/window;
      Nf{cnt1} = Nf{cnt1}/window;
      
      Xf{cnt1} = Xf{cnt1}(1:window/2+1,:);
      Yf{cnt1} = Yf{cnt1}(1:window/2+1,:);
      Sf{cnt1} = Sf{cnt1}(1:window/2+1,:);
      Nf{cnt1} = Nf{cnt1}(1:window/2+1,:);
      
      Xf{cnt1}(2:end-1,:) = 2*Xf{cnt1}(2:end-1,:);
      Yf{cnt1}(2:end-1,:) = 2*Yf{cnt1}(2:end-1,:);
      Sf{cnt1}(2:end-1,:) = 2*Sf{cnt1}(2:end-1,:);
      Nf{cnt1}(2:end-1,:) = 2*Nf{cnt1}(2:end-1,:);

      xy = sum(Xf{cnt1} .* conj(Yf{cnt1}), 2)/N;
      xx = sum(Xf{cnt1} .* conj(Xf{cnt1}), 2)/N;
      yy = sum(Yf{cnt1} .* conj(Yf{cnt1}), 2)/N;

      ss = sum(Sf{cnt1} .* conj(Sf{cnt1}), 2)/N;
      nn = sum(Nf{cnt1} .* conj(Nf{cnt1}), 2)/N;

      coh_temp(cnt1,:) = sqrt(abs(xy).^2./(xx.*yy));
      snr_temp(cnt1,:) = ss./nn;

      Pow_temp1(cnt1,:) = sum(abs(Xf{cnt1}),2)/N;
      Pow_temp2(cnt1,:) = sum(abs(Yf{cnt1}),2)/N; 
      S_temp(cnt1,:) = sum(abs(Sf{cnt1}),2)/N;
      N_temp(cnt1,:) = sum(abs(Nf{cnt1}),2)/N;  

end

Coh  = sqrt(abs(mean(xy,1)).^2./(mean(xx,1).*mean(yy,1)));
SNR = mean(ss,1)./mean(nn,1);

Coh = Coh(1:window/2+1);
SNR = SNR(1:window/2+1);

pow(1,:) = mean(Pow_temp1,1);
pow(2,:) = mean(Pow_temp2,1);
pow = pow(:,1:window/2+1);

S(1,:) = mean(S_temp,1);
Ns(2,:) = mean(N_temp,1);

S = S(:,1:window/2+1);
Ns = Ns(:,1:window/2+1);


fout = params.fsample*(0:(window/2))/window;

if isfield(params,'freq')
    [~, i] = min(abs(fout - mean(params.freq)));
    SNRf = SNR(i);
    Cohf = Coh(i);
else
    [val id] = max(S);
    SNRf = SNR(id);
    Cohf = Coh(id);    
end

end