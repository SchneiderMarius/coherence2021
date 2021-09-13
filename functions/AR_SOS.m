function [fout,SNR,SNRf,S,Ns] = AR_CohSOS(Signal,Noise,foi,params,window,ovl)

if nargin<6
    window = 1000;
    ovl = window/2;
elseif nargin<7
    ovl = window/2;    
end

coh_temp = [];  

if iscell(Signal)
    for cnt1 = 1 : length(Signal)
      Signal{cnt1} = Signal{cnt1}(1:end,1);
      Noise{cnt1}  = Noise{cnt1}(1:end,1);


    %   xt{cnt1} = reshape(x{cnt1}(1:L*N), L, N);
    %   yt{cnt1} = reshape(y{cnt1}(1:L*N), L, N); 
    %   St{cnt1} = reshape(Signal{cnt1}(1:L*N), L, N);
    %   Nt{cnt1} = reshape(Noise{cnt1}(1:L*N), L, N);  

      z        = buffer(1:numel(Signal{cnt1}),window,ovl);
      St{cnt1} = Signal{cnt1}(z(:,all(z)));  

      z        = buffer(1:numel(Noise{cnt1}),window,ovl);
      Nt{cnt1} = Noise{cnt1}(z(:,all(z)));    

      N = size(Nt{cnt1},2);

      % transform to frequency domain
      Sf{cnt1} = fft(St{cnt1},window,1);
      Nf{cnt1} = fft(Nt{cnt1},window,1);    

      Sf{cnt1} = Sf{cnt1}/window;
      Nf{cnt1} = Nf{cnt1}/window;
      
      Sf{cnt1} = Sf{cnt1}(1:window/2+1,:);
      Nf{cnt1} = Nf{cnt1}(1:window/2+1,:);
      
      Sf{cnt1}(2:end-1,:) = 2*Sf{cnt1}(2:end-1,:);
      Nf{cnt1}(2:end-1,:) = 2*Nf{cnt1}(2:end-1,:);
          

      % estimate expectations by taking the average over N blocks
      ss(cnt1,:) = sum(Sf{cnt1} .* conj(Sf{cnt1}), 2)/N;
      nn(cnt1,:) = sum(Nf{cnt1} .* conj(Nf{cnt1}), 2)/N;
        
      S_temp(cnt1,:)    = sum(abs(Sf{cnt1}),2)/N;
      N_temp(cnt1,:)    = sum(abs(Nf{cnt1}),2)/N;  

    end
elseif ~iscell(Signal)
      cnt1 = 1;  
      Signal2{cnt1} = Signal(101:end,1);
      Noise2{cnt1}  = Noise(101:end,1);

      z        = buffer(1:numel(x{cnt1}),window,ovl);
      xt{cnt1} = x{cnt1}(z(:,all(z)));   

      z        = buffer(1:numel(Noise2{cnt1}),window,ovl);
      Nt{cnt1} = Noise2{cnt1}(z(:,all(z)));    

      N = size(Nt{cnt1},2);

      % transform to frequency domain
      Sf{cnt1} = fft(St{cnt1},window,1);
      Nf{cnt1} = fft(Nt{cnt1},window,1);    

      Sf{cnt1} = Sf{cnt1}/window;
      Nf{cnt1} = Nf{cnt1}/window;
      
      Sf{cnt1} = Sf{cnt1}(1:window/2+1,:);
      Nf{cnt1} = Nf{cnt1}(1:window/2+1,:);
      
      Sf{cnt1}(2:end-1,:) = 2*Sf{cnt1}(2:end-1,:);
      Nf{cnt1}(2:end-1,:) = 2*Nf{cnt1}(2:end-1,:);

      % estimate expectations by taking the average over N blocks
      ss = sum(Sf{cnt1} .* conj(Sf{cnt1}), 2)/N;
      nn = sum(Nf{cnt1} .* conj(Nf{cnt1}), 2)/N;

      snr_temp(cnt1,:) = ss./nn;

      S_temp(cnt1,:) = sum(abs(Sf{cnt1}),2)/N;
      N_temp(cnt1,:) = sum(abs(Nf{cnt1}),2)/N;  

end

% coh_temp(cnt1,:)  = sqrt(abs(mean(xy,2)).^2./(mean(xx,2).*mean(yy,2)));
% 
% figure
% plot(fout,coh_temp(cnt1,:))
SNR = mean(ss,1)./mean(nn,1);

% Coh = mean(coh_temp,1);
SNR = SNR(1:window/2+1);

S(1,:) = mean(S_temp,1);
Ns(1,:) = mean(N_temp,1);

S = S(:,1:window/2+1);
Ns = Ns(:,1:window/2+1);


fout = params.fsample*(0:(window/2))/window;

if isfield(params,'freq')
    [~, i] = min(abs(fout - mean(params.freq)));
    SNRf = SNR(i);
else
    [val id] = max(S);
    SNRf = SNR(id);
end

end