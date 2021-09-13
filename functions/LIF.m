function [V,spt,Frate] = LIF(I,ge,params,N)

if nargin<2 || isempty(ge)
    ge = zeros(length(I),N);
end

if nargin<3 || isempty(params)
   params.tau   = 20;
   params.Vrest = -60;
   params.Vreset = -50;
   params.R     = 10;
   params.thres = -30;
   params.Ee    = 60;
   params.amp = 40;
end

params.dt    = 1;
% params.tstop = 200;
params.tstop = length(I)*params.dt;

spt = zeros(params.tstop/params.dt,1);

V(1,1:N) = params.Vrest;
spt      = sparse(zeros(length(I),N));

if size(I,2)<N
    I(:,1:N) = repmat(I(:,1),1,N);
end


for cnt = 1 : params.tstop/params.dt-1
    
    V(cnt+1,:) = params.dt/params.tau*(-V(cnt,:)+params.Vrest + params.R*I(cnt,:) + ge(cnt)) + V(cnt,:);
    
    if any(V(cnt+1,:)>params.thres)
       id = V(cnt+1,:)>params.thres;
       V(cnt,id) = params.amp;
       V(cnt+1,id) = params.Vreset;
       spt(cnt+1,id) = 1;
    end
end

Frate = sum(spt)/(params.tstop/1000);
% if plot_it==1
%     figure
%     subplot(2,1,1)
%     hold all
%     plot(params.dt:params.dt:params.tstop,V);
%     title('voltage trace')
%     subplot(2,1,2)
%     hold all
%     for cnt1 = 1 : N
%        plot(params.dt:params.dt:params.tstop,spt(:,cnt1)*cnt1,'k.');
%     end
%     ylim([0.5 N+0.5])
%     
%     title('input current')
% %     plot(params.dt:params.dt:params.tstop,I);
% %     subplot(3,1,3)
% %     plot(params.dt:params.dt:params.tstop,ge);
% %     title('synaptic current')
% end
end