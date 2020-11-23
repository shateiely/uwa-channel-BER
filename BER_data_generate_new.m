% 用来产生调制和信道的关系的数据
% 信道使用随机信道
% 调制方式为 BPSK QPSK 8PSK
% 通信体制为单载波频域均衡
% 接收端首先假设信道已知

clc
close all
clear


Len_block = 1024; % 包含CP和数据的总长度


chan_order=50;
N_cp=256;

PN=randi([0,1], 1, 512);    % 1x512，取值0-1的伪随机矩阵-->长512的伪随机码
temp=[PN(1,1:N_cp)   ]; % 取PN的前一半(256)
CP=temp+1j*temp;    % 1->1+j,0->0+0j
N_cp =length(CP);   % cp的长度
N_data = Len_block-N_cp; % 数据长度=CP和数据的总长度-CP长度 #of symbols in a fde block



NumTrials = 1; % #of frames in each evaluation

nfft1=N_data;          % 原始数据块长度
nfft_PN=nfft1+N_cp;   % 采用PN处理的fft长度


%
SNR = 12:2:16; % 三档信噪比12，14，16
SER = zeros(4,length(SNR)); % zf, mmse, dfe, ib-HD锛宨b-sd 横4竖3
nblock =10; % fde blocks in each frame, >1考虑IBI; <=1无IBI
Len_chan=30;
M=3; % 调制阶数
Num_chan=10000;
SER_MMSE=zeros(Num_chan,M,length(SNR)); % 10000x3x3
SER_MMSE_mean=zeros(Num_chan,M,length(SNR));

N_SNR=length(SNR);

rng('default') % 初始化随机数生成器，使结果具备可重复性。
N_data = (Len_block-N_cp)/M; % 数据长度768/M3=256
x =zeros(N_data,M); % 256行3列表格
xtx_pload = zeros(N_data,M); % 256行3列表格
for mod =1:M
    M_mod=2^mod; % 2 4 8
    x(:,mod) = randi([0,M_mod-1], N_data,1); % 在第（阶数）列赋给256x1，取值0-1、3、7的伪随机列，
    xtx_pload(:,mod) =pskmod(x(:,mod),M_mod) ; % 第一列BPSK，第二列QPSK，第三列8PSK
end
xtx_pload = reshape(xtx_pload,[],1); % 256x3的xtx_pload重构为768x1的矩阵
Pilot = [CP.';xtx_pload;CP.'];  % 加导频256+768+256
Len_x = length(Pilot); % 1280
train_data=zeros(Num_chan*M*length(SNR)*NumTrials,Len_x*2+1); % 10000*3*3*1，2561
BER_target=zeros(Num_chan*M*length(SNR)*NumTrials,1); % 10000*3*3*1

for index_chan=1:Num_chan % 1到10000
    
    rng(index_chan)
    taps=randi([1,Len_chan],1,1); %1到30随机一个数
    chan=Generat_Channel(taps,Len_chan,index_chan); % 信道：1x30的复数（随机的taps，30，for1—10000）
    chan_add=zeros(1,chan_order); % 1x50的表格
    chan_add(1:length(chan))=chan; % chan的1x30复数赋给chan_add的1x50
    for mod=1:M
        M_mod=2^mod; % 三种调制方式 BPSK QPSK 8PSK
        for n = 1:length(SNR) % 三档信噪比 14 16 18
            for m = 1:NumTrials % quasi-static channel modeling 时不变信道建模
                
                x = randi([0,M_mod-1], N_data, nblock); % 256x10 0到1、3、7 N_data = 256 nblock =10
                xtx_pload =pskmod(x,M_mod) ; % 2560的数据调制成xPSK
                %         xtx_pload=ones(128,100);
                xtx_pload_Power = sum(abs(xtx_pload).^2,1)/size(xtx_pload,1); % 发射功率
                
                PN_add=[]; % PN:伪随机码
                for i=1:nblock
                    PN_add=[PN_add CP.']; % 前一半PN+j前一半PN
                end
                xtx = [xtx_pload;PN_add];  %% 添加PN
                xtx_temp= reshape(xtx,[],1); % 重构为?x1的矩阵
                xtx_all= [Pilot;xtx_temp];  % 加导频
                index=(index_chan-1)*M*N_SNR*NumTrials+(mod-1)*M*NumTrials+(n-1)*NumTrials+m;
                % 选择train_data中对应信噪比的位置; index_chan=1到1w，M=3，N_SNR=SNR长度
                train_data(index,1)=mod;
                fadesig = filter(chan,1, xtx_all); % Effect of channel, quasi-static
                
                %%
                xtxPower = sum(abs(xtx(:)).^2)/length(xtx(:)); % 信号能量
                xtxPower_dB = 10*log10(xtxPower); % 转换为dB
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
                
                h_CE=CE_IPNLMS(rnoise(1:N_cp),CP,chan_order); % 信道估计
                train_data(index,2:1+length(rnoise_input))=real(rnoise_input); % 估计的信道 实部
                train_data(index,1+length(rnoise_input)+1:1+length(rnoise_input)*2)=imag(rnoise_input); % 虚部
                
                
                rnoise = rnoise(Len_x+1:end, :); % 移除PN
                rnoise_CP = reshape(rnoise, [], nblock);  %%
                
                % H=fft(chan, nfft_PN, 2); % 已知信道
                H=fft(h_CE, nfft_PN, 2); % 未知信道
                
                for i_block=1:nblock
                    
                    rnoise_inter=rnoise_CP(:,i_block);  %% 取一个数据块
                    fdein = fft(rnoise_inter,nfft_PN); % 转换为FDE的频域进行处理
                    
                    H_mmse = (H')./((abs(H).^2+ones(size(H))/(xtxPower/noisePower)).');% mmse，发射信号与接收器的噪声功率比
                    fdeout_mmse = fdein.*H_mmse;
                    xrx_mmse = ifft(fdeout_mmse,nfft_PN); % 转换回时域
                    % z_mmse=demodulate(modem.qamdemod(M_mod),xrx_mmse);
                    z_mmse = pskdemod(xrx_mmse,M_mod); % psk解码
                    SER_MMSE_temp(i_block)=BER_Cacula(z_mmse(1:N_data),x(:,i_block));
                    % SER计算
                    
                end
                BER_target(index)=mean(SER_MMSE_temp); % BER=SER均方
            end
        end
        %         SER_MMSE_mean(index_chan,mod,n) = SER_MMSE(index_chan,mod,n) /NumTrials/nblock;
        %         %             SER_FDO_MMSE_mean = SER_FDO_MMSE/numTrials/nblock;
        %
    end
    
    % close(h_waitbar)
    
    %     SER_FDO_MMSE_CE_mean= SER_FDO_MMSE_CE/numTrials/nblock;
    %     SER_MMSE_CE_mean= SER_MMSE_CE/numTrials/nblock;
end
%%

% save SER_MMSE_mean SER_MMSE_mean
save('Data_set_all', 'train_data', 'BER_target');
