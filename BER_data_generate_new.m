% 用来产生调制和信道的关系的数据
% 信道使用随机信道
% 调制方式为 BPSK QPSK 8PSK
% 通信体制为单载波频域均衡
% 接收端首先假设信道已知

clc
close all
clear


Len_block = 1024; % 包含CP和数据的总长度
chan_order=30;
N_cp=256;
N_data = Len_block-N_cp;

PN=2*randi([0,1], 1, 512)-1;    % 1x512，取值0-1的伪随机矩阵-->长512的伪随机码
CP=PN(1,1:N_cp); % 取部分伪随机序列
N_cp =length(CP);   % cp的长度
NumTrials = 1; % #of frames in each evaluation

nfft_PN=Len_block;   % 采用PN处理的fft长度

SNR = 12:2:16; % 三档信噪比12，14，16
nblock =10; % fde blocks in each frame, >1考虑IBI; <=1无IBI
Len_chan=30;    %信道时延扩展 
M=3; % 调制阶数
Num_chan=10000;

N_SNR=length(SNR);

rng('default') % 初始化随机数生成器，使结果具备可重复性。
N_pilot = (Len_block-N_cp)/M; % 每种调制256个导频符号
x =zeros(N_pilot,M); % 256行3
xtx_pload = zeros(N_pilot,M); % 256行3列表格
for mod =1:M
    M_mod=2^mod; % 2 4 8
    x(:,mod) = randi([0,M_mod-1], N_pilot,1);
    xtx_pload(:,mod) =pskmod(x(:,mod),M_mod);
end
xtx_pload = reshape(xtx_pload,[],1); % 256x3重构为768x1
Pilot = [CP.';xtx_pload;CP.'];  % 加循环前缀256+768+256
Len_x = length(Pilot); % 1280

PN_add=zeros(N_cp,nblock);
for i=1:nblock
    PN_add(:,i)=CP.'; % CP:256x1,nblock个CP
end
train_data=zeros(Num_chan*M*length(SNR)*NumTrials,Len_x*2+1); % 10000*3*3*1，2561
BER_target=zeros(Num_chan*M*length(SNR)*NumTrials,1); % 10000*3*3*1

for index_chan=1:Num_chan
    
    rng(index_chan)
    taps=randi([1,Len_chan],1,1);
    chan=Generat_Channel(taps,Len_chan,index_chan); % 信道：1x30的复数（随机的taps，30，for1—10000）
    for n = 1:length(SNR) % 三档信噪比 14 16 18
        for mod=1:M
            M_mod=2^mod; % 三种调制方式 BPSK QPSK 8PSK
            for m = 1:NumTrials % quasi-static channel modeling 时不变信道建模
                
                x = randi([0,M_mod-1], N_data, nblock); % 256x10 0到1、3、7 N_data = 256 nblock =10
                xtx_pload =pskmod(x,M_mod) ; % 2560的数据调制成xPSK
                %         xtx_pload=ones(128,100);
                xtx_pload_Power = sum(abs(xtx_pload).^2,1)/size(xtx_pload,1); % 发射功率
                
                xtx = [xtx_pload;PN_add];  %% 添加PN
                xtx_temp= reshape(xtx,[],1); % 重构为?x1的矩阵
                xtx_all= [Pilot;xtx_temp];  % 加导频
                index=(index_chan-1)*M*N_SNR*NumTrials+(n-1)*M*NumTrials+(mod-1)*NumTrials+m;
                % 选择train_data中对应信噪比的位置; index_chan=1到1w，M=3，N_SNR=SNR长度
                train_data(index,1)=mod;
                fadesig = filter(chan,1, xtx_all); % Effect of channel, quasi-static
                
                %%
                xtxPower = sum(abs(xtx(:)).^2)/length(xtx(:)); % 信号能量
                xtxPower_dB = 10*log10(xtxPower);
                fadesigPower = sum(abs(fadesig(:)).^2)/length(fadesig(:)); % 衰落信号能量
                fadesigPower_dB = 10*log10(fadesigPower); % 转换为dB
                %         noisePower_dB = fadesigPower_dB-SNR(n); % compute noise power, in dB
                noisePower_dB = xtxPower_dB-SNR(n);       % 噪声能量
                noisePower = 10^(noisePower_dB/10); % 转换为线性
                rng(NumTrials) % 使用NumTrials为随机数生成器提供种子，以使 rand、randi 和 randn 生成可预测的数字序列。
                noise=sqrt(noisePower)*sqrt(1/2)*(randn(size(fadesig))+1i*randn(size(fadesig)));
                % rnoise = awgn(fadesig,SNR(n),'measured'); % Additive receiver noise, receiver snr
                % rnoise = awgn(fadesig, SNR(n), xtxPower_dB); % additive noise with the assumption
                rnoise  =  fadesig+ noise;
                rnoise_input = rnoise(1:Len_x);
                
                h_CE=CE_IPNLMS(rnoise(Len_x-N_cp+1:Len_x),CP,chan_order); % 一个时刻的信道估计
                train_data(index,2:1+length(rnoise_input))=real(rnoise_input); % 估计的信道 实部
                train_data(index,1+length(rnoise_input)+1:1+length(rnoise_input)*2)=imag(rnoise_input); % 虚部
                
                
                rnoise = rnoise(Len_x+1:end, :); % 移除pilot
                rnoise_CP = reshape(rnoise, [], nblock);  %%
                
                % H=fft(chan, nfft_PN, 2); % 已知信道
                
                SER_MMSE_temp =zeros(1,nblock);
                for i_block=1:nblock
                    H=fft(h_CE, nfft_PN, 2); % 未知信道
                    rnoise_inter=rnoise_CP(:,i_block);  %% 取一个数据块
                    fdein = fft(rnoise_inter,nfft_PN); % 转换到频域进行处理
                    
                    H_mmse = (H')./((abs(H).^2+ones(size(H))/(xtxPower/noisePower)).'); 频域均衡
                    fdeout_mmse = fdein.*H_mmse;
                    xrx_mmse = ifft(fdeout_mmse,nfft_PN); % 转换回时域
                    z_mmse = pskdemod(xrx_mmse,M_mod); % psk解码
                    SER_MMSE_temp(i_block)=BER_Cacula(z_mmse(1:N_data),x(:,i_block));
                    % SER计算
                    h_CE=CE_IPNLMS(rnoise_inter(end-N_cp+1:end),CP,chan_order); % 信道估计
                end
                BER_target(index)=mean(SER_MMSE_temp); % BER=SER均方
            end
        end
    end

end

save('Data_set_all', 'train_data', 'BER_target');
