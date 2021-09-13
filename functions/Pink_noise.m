function pink = Pink_Noise(len,num,a)

rng('shuffle')

pink=[];
for cnt = 1 : num
    x = (randn(1,len));
    X = fft(x);
    numunique = len/2+1;
    k = 1 : numunique;
    X = X(:,1:numunique);
    X = X./(k.^(a));
    X = [X conj(X(end-1:-1:2))];
    pink(:,cnt) = real(ifft(X));
end

% Pow = abs(fft(x)./st);
% Pow = Pow(1:st/2+1);
% Pow = 2*Pow(1:end);
% 
% freq = (fs)*(0:(st/2))/st;
% 
% plot(freq,Pow)
end