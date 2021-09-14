function [parsf,outpf] = AR2model(freq,fsample,eigenvalue,SampleSize)

% find coefficients for AR-2 model oscillating at freq with specific
% eigenvalue

if nargin<4
   SampleSize = 2000000;
end

cnt = 0;

outp = []; pars = []; ang2 = [];
for k = 1:SampleSize
    coeff = -2+4*rand(1,2);
    rt = roots([-coeff 1]);
    z = 1./rt;
    if ~isreal(rt)        
        mag = abs(z);
        if mag(1)<1 && mag(2)<1
            cnt = cnt + 1;
            ang = atan2(imag(z(1)), real(z(1)));
            outp(cnt,:) = [mag(1) abs(ang)];
            pars(cnt,:) = coeff;
            ang2(cnt) = acos((pars(cnt,1).*(pars(cnt,2)-1))./(4.*pars(cnt,2)))./(2*pi);
        end
    end
end
   

v1 = 2*pi*freq(1)./fsample; 
v2 = 2*pi*freq(2)./fsample;

indx = find(outp(:,2)>v1 & outp(:,2)<v2);

outpf = outp(indx,:);
parsf = pars(indx,:);

if nargin==3 && ~isempty(eigenvalue)
    dist    = abs(outpf(:,1)-eigenvalue);  
    minDist = min(dist);
	ind     = find(dist == minDist);
    parsf = parsf(ind,:);
    outpf = outpf(ind,:);
end

end

