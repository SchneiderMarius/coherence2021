function [fout,Coh,Cohf,pow] = AR_Coh(LFP,foi,params,window,ovl)

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

    %   xt{cnt1} = reshape(x{cnt1}(1:L*N), L, N);
    %   yt{cnt1} = reshape(y{cnt1}(1:L*N), L, N); 
    %   St{cnt1} = reshape(Signal{cnt1}(1:L*N), L, N);
    %   Nt{cnt1} = reshape(Noise{cnt1}(1:L*N), L, N);  

      z        = buffer(1:numel(x{cnt1}),window,ovl);
      xt{cnt1} = x{cnt1}(z(:,all(z)));   

      z        = buffer(1:numel(y{cnt1}),window,ovl);
      yt{cnt1} = y{cnt1}(z(:,all(z)));     

      N = size(yt{cnt1},2);

      % transform to frequency domain
      Xf{cnt1} = fft(xt{cnt1},window,1);
      Yf{cnt1} = fft(yt{cnt1},window,1);   
   

      Xf{cnt1} = Xf{cnt1}/window;
      Yf{cnt1} = Yf{cnt1}/window;

      
      Xf{cnt1} = Xf{cnt1}(1:window/2+1,:);
      Yf{cnt1} = Yf{cnt1}(1:window/2+1,:);

      Xf{cnt1}(2:end-1,:) = 2*Xf{cnt1}(2:end-1,:);
      Yf{cnt1}(2:end-1,:) = 2*Yf{cnt1}(2:end-1,:);
         

      % estimate expectations by taking the average over N blocks
      xy(cnt1,:) = sum(Xf{cnt1} .* conj(Yf{cnt1}), 2)/N;
      xx(cnt1,:) = sum(Xf{cnt1} .* conj(Xf{cnt1}), 2)/N;
      yy(cnt1,:) = sum(Yf{cnt1} .* conj(Yf{cnt1}), 2)/N;

        
      Pow_temp1(cnt1,:) = sum(abs(Xf{cnt1}),2)/N;
      Pow_temp2(cnt1,:) = sum(abs(Yf{cnt1}),2)/N;

    end
elseif ~iscell(LFP)
      cnt1 = 1;  
      x{cnt1}    = LFP(101:end,1);
      y{cnt1}    = LFP(101:end,2);


      z        = buffer(1:numel(x{cnt1}),window,ovl);
      xt{cnt1} = x{cnt1}(z(:,all(z)));   

      z        = buffer(1:numel(y{cnt1}),window,ovl);
      yt{cnt1} = y{cnt1}(z(:,all(z)));     


      N = size(yt{cnt1},2);

      % transform to frequency domain
      Xf{cnt1} = fft(xt{cnt1},window,1);
      Yf{cnt1} = fft(yt{cnt1},window,1);   

 

      Xf{cnt1} = Xf{cnt1}/window;
      Yf{cnt1} = Yf{cnt1}/window;

      
      Xf{cnt1} = Xf{cnt1}(1:window/2+1,:);
      Yf{cnt1} = Yf{cnt1}(1:window/2+1,:);

      Xf{cnt1}(2:end-1,:) = 2*Xf{cnt1}(2:end-1,:);
      Yf{cnt1}(2:end-1,:) = 2*Yf{cnt1}(2:end-1,:);


      % estimate expectations by taking the average over N blocks
      xy = sum(Xf{cnt1} .* conj(Yf{cnt1}), 2)/N;
      xx = sum(Xf{cnt1} .* conj(Xf{cnt1}), 2)/N;
      yy = sum(Yf{cnt1} .* conj(Yf{cnt1}), 2)/N;

      coh_temp(cnt1,:) = sqrt(abs(xy).^2./(xx.*yy));

      Pow_temp1(cnt1,:) = sum(abs(Xf{cnt1}),2)/N;
      Pow_temp2(cnt1,:) = sum(abs(Yf{cnt1}),2)/N; 
end

% coh_temp(cnt1,:)  = sqrt(abs(mean(xy,2)).^2./(mean(xx,2).*mean(yy,2)));
% 
% figure
% plot(fout,coh_temp(cnt1,:))
Coh  = sqrt(abs(mean(xy,1)).^2./(mean(xx,1).*mean(yy,1)));

% Coh = mean(coh_temp,1);
Coh = Coh(1:window/2+1);

pow(1,:) = mean(Pow_temp1,1);
pow(2,:) = mean(Pow_temp2,1);
pow = pow(:,1:window/2+1);

fout = params.fsample*(0:(window/2))/window;

if isfield(params,'freq')
    [~, i] = min(abs(fout - mean(params.freq)));
    Cohf = Coh(i);
else
    [val id] = max(S);
    Cohf = Coh(id);    
end

end