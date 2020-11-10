% 产生随机信道，
% Len_chan： 信道长度
% taps: 信道多径个数
% index_chan: 信道索引
function chan=Generat_Channel(taps,Len_chan,index_chan)
rng('default')
rng(index_chan)
a=randperm(Len_chan);
delay_taps=a(1:taps);
pow_prof = (1/taps) * (ones(1,taps));
chan_coef = sqrt(pow_prof).*(sqrt(1/2) * (randn(1,taps)+1i*randn(1,taps)));
chan=zeros(1,Len_chan);
chan(delay_taps)=chan_coef;