function [h]=CE_IPNLMS(rec_sig,training,chan_order)
%  % [h]=CE_IPNLMS(rec_sig,training,W,chan_order)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 输入参数：
% rec_sig:  经过信道和噪声的接收信号
% training: 训练序列
% chan_order:  信道长度
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 输出参数：
% h: 信道
%
%   Author(s): 何成兵

%   Copyright 2016-2035 西北工业大学.
%   $Revision: 1.0 $  $Date: 2016/05/13   $
[Nuser,~]=size(training);
W  = zeros(chan_order*Nuser,1);
% W  = ones(chan_order,1)/sqrt(126);
%%% IPNLMS的参数
alpha =0;kezi =0.01; deta =0.01;
stepa = 0.4;  %%%注意和LMS的miu不一样，miu要尽量小
trainlen=length(training) ;


AA = [zeros(chan_order-1,Nuser); training.'; zeros(chan_order,Nuser)];
weif=zeros(1,chan_order);
Kl=zeros(1 ,chan_order*Nuser);
for n=1:trainlen
    
    U = AA(chan_order+n-1:-1:n,:);
    U2=reshape(U.',[],1);
    y = W' * U2;
    d = rec_sig(n);
    e = d - y;
    Norm_1 = sum(abs(W));
    
    for l =1 :chan_order*Nuser
        Kl(l) = (1-alpha)/(2*(chan_order))+(1+alpha)*abs(W(l))/(2*Norm_1+kezi);
    end
    K = diag(Kl);
    %     W1= W + stepa*K*e(n)'*U/(U'*K*U+deta);
    W= W + stepa*e'*Kl.'.*U2/(U2'.*Kl*U2+deta);
    %     plot(abs(W'))
    % pause(0.1)
end

h=W';



