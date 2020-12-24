% 用来产生调制和信道的关系的数据
% 信道使用随机信道
% 调制方式为 BPSK QPSK 8PSK
% 通信体制为单载波频域均衡
% 接收端首先假设信道已知

clc
close all
clear

Tc=0.25*1e-3; % 符号间隔
                
                
Len_block = 1024; % 包含CP和数据的总长度


chan_order=120;
N_cp=512;

PN=randi([0,1], 1, 512); % 1x512，取值0-1的伪随机矩阵-->长512的伪随机码
temp=[PN(1,1:N_cp)   ]; % 取PN的前一半(256)
CP=temp+1j*temp; % 1->1+1j,0->0+0j
N_cp =length(CP); % cp的长度
N_data = Len_block-N_cp; % 数据长度=CP和数据的总长度-CP长度

NumTrials = 1; % #of frames in each evaluation

nfft1=N_data;          % 原始数据块长度
nfft_PN=nfft1+N_cp;   % 采用PN处理的fft长度

%
SNR = 12:2:16; % 三档信噪比12，14，16
SER = zeros(4,length(SNR)); % zf, mmse, dfe, ib-HD锛宨b-sd
nblock =10; % >1考虑IBI; <=1无IBI
Len_chan=100; % 信道多径个数 1:30
M=3; % 调制阶数
Num_chan=156;
SER_MMSE=zeros(Num_chan,M,length(SNR));
SER_MMSE_mean=zeros(Num_chan,M,length(SNR));

N_SNR=length(SNR);

rng('default') % 初始化随机数生成器，使结果具备可重复性。
N_data_pilot = round((Len_block-N_cp)/M); % 导频长度768/M3=256
x_pilot =zeros(N_data_pilot,M);
xtx_pload = zeros(N_data_pilot,M);

for mod =1:M
    M_mod=2^mod;
    x_pilot(:,mod) = randi([0,M_mod-1], N_data_pilot,1);
    xtx_pload(:,mod) =pskmod(x_pilot(:,mod),M_mod) ;
end
xtx_pload = reshape(xtx_pload,[],1);
Pilot = [CP.';xtx_pload;CP.'];  %% add PN
Len_x = length(Pilot);
train_data=zeros(Num_chan*M*length(SNR)*NumTrials,Len_x*2+1);
BER_target=zeros(Num_chan*M*length(SNR)*NumTrials,1);

for index_chan=1:Num_chan
    
    rng(index_chan)
    
    for mod=1:M
        M_mod=2^mod; % 三种调制方式
        for n = 1:length(SNR) % 三档信噪比 14 16 18
            for m = 1:NumTrials % quasi-static channel modeling
                
                x = randi([0,M_mod-1], N_data, nblock);
                xtx_pload =pskmod(x,M_mod) ;
                %         xtx_pload=ones(128,100);
                xtx_pload_Power = sum(abs(xtx_pload).^2,1)/size(xtx_pload,1); % check transmitted signal power along each tx antenna
                
                PN_add=[];
                for i=1:nblock
                    PN_add=[PN_add CP.']; % CP:256x1,nblock个CP
                end
                xtx = [xtx_pload;PN_add];  %% add PN
                xtx_temp= reshape(xtx,[],1);
                xtx_all= [Pilot;xtx_temp];
                index=(index_chan-1)*M*N_SNR*NumTrials+(mod-1)*M*NumTrials+(n-1)*NumTrials+m;
                %             train_data(index,1)=n;
                train_data(index,1)=mod;
                % fadesig = filter(chan,1, xtx_all); % Effect of channel, quasi-static
% load(('KAU1_1.mat'), 'h')
% h=h/abs(max(h(:))); % 归一化
% resample_multiple=4; % 抽样倍数
% resample_path=1;
% for h_path=1:size(h,2)
%     if rem(h_path,resample_multiple)==0
%         h_resample(:,resample_path)=h(:,h_path);
%         resample_path=resample_path+1;
%     end
% end
% Channel_temp=h_resample(:,1:Len_chan).';
% sampling_temp=size(Channel_temp,2);
% interp_multiple=ceil(length(xtx_all)/sampling_temp); % 插值倍数
% for nn=1:Len_chan
%     channel_interp(nn,:) = interp(Channel_temp(nn,:),interp_multiple); %对每列进行multiple倍插值
% end
% Channel_Real_temp=channel_interp(:,1:length(xtx_all)); % 截取与发射信号等长的信道
% Channel_Real=Channel_Real_temp/abs(max(Channel_Real_temp(:)));

Channel_Real=Choose_channel(Len_chan,index_chan,xtx_all);
                [N_path,N_sampling]=size(Channel_Real);
                Rxtemp=zeros(N_path,N_path+length(xtx_all)-1);
                for itao=1:N_path
                    Rxtemp(itao,itao:length(xtx_all)+itao-1) =  Channel_Real(itao,:).*xtx_all.';
                end
                fadesig_TV_temp= sum(Rxtemp,1);
                fadesig = fadesig_TV_temp(1:length(xtx_all)).';
                
                %%
                xtxPower = sum(abs(xtx(:)).^2)/length(xtx(:)); % check transmitted signal power
                xtxPower_dB = 10*log10(xtxPower); % convert to dB
                fadesigPower = sum(abs(fadesig(:)).^2)/length(fadesig(:)); % check fading signal power
                fadesigPower_dB = 10*log10(fadesigPower); % convert to dB
                %         noisePower_dB = fadesigPower_dB-SNR(n); % compute noise power, in dB
                noisePower_dB = xtxPower_dB-SNR(n);       % compute noise power, with
                noisePower = 10^(noisePower_dB/10); % convert backto linear scale
                rng(NumTrials)
                noise=sqrt(noisePower)*sqrt(1/2)*(randn(size(fadesig))+1i*randn(size(fadesig)));
                %         rnoise = awgn(fadesig,SNR(n),'measured'); % Additive receiver noise, receiver snr
                %         rnoise = awgn(fadesig, SNR(n), xtxPower_dB); % additive noise with the assumption
                
                %rnoise  =  fadesig;
                rnoise  =  fadesig+ noise;
                rnoise_input = rnoise(1:Len_x);
                
                h_CE=CE_IPNLMS(rnoise(1:N_cp),CP,chan_order);
                
                train_data(index,2:1+length(rnoise_input))=real(rnoise_input); % 前一半是实部
                train_data(index,1+length(rnoise_input)+1:1+length(rnoise_input)*2)=imag(rnoise_input); % 后一半是虚部
                
                
                rnoise = rnoise(Len_x+1:end, :); % remove PN
                rnoise_CP = reshape(rnoise, [], nblock);  %%
                % H=fft(chan, nfft_PN, 2); % 已知信道
                H=fft(h_CE, nfft_PN, 2); % 未知信道
                
                for i_block=1:nblock
                    
                    rnoise_inter=rnoise_CP(:,i_block);  %% 鍙栦竴涓寘
                    fdein = fft(rnoise_inter,nfft_PN); % convert into frequency domain for FDE
                    
                    H_mmse = (H')./((abs(H).^2+ones(size(H))/(xtxPower/noisePower)).');% mmse with transmitted signal to receiver noise power ratio
                    fdeout_mmse = fdein.*H_mmse;
                    xrx_mmse = ifft(fdeout_mmse,nfft_PN); % convert back to time domain for detection
                    %                 z_mmse=demodulate(modem.qamdemod(M_mod),xrx_mmse);
                    z_mmse = pskdemod(xrx_mmse,M_mod);
                    SER_MMSE_temp(i_block)=BER_Cacula(z_mmse(1:N_data),x(:,i_block));
                    %                     SER_MMSE(index_chan,mod,n) =SER_MMSE(index_chan,mod,n) +SER_MMSE_temp;
                    
                end
                BER_target(index)=mean(SER_MMSE_temp);
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
save('Data_set_Real', 'train_data', 'BER_target');
