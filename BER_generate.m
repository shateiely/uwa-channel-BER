% 用来产生调制和信道的关系的数据
% 信道首先使用随机信道
% 调制方式为 BPSK QPSK 8PSK
% 通信体制为单载波频域均衡
% 接收端首先假设信道已知

clc
close all
clear
rng('default')


Len_block = 512; % 包含CP和数据的总长度



chan_order=50;
N_cp=256;

PN=randi([0,1], 1, 512);
temp=[PN(1,1:N_cp)   ];
CP=temp+1j*temp;
N_cp =length(CP);
N_data = Len_block-N_cp; % #of symbols in a fde block



numTrials = 100; % #of frames in each evaluation

nfft1=N_data;          % 原始数据块长度
nfft_PN=nfft1+N_cp;   % 采用PN处理的fft长度


%
SNR = 8:2:12;
SER = zeros(4,length(SNR)); % zf, mmse, dfe, ib-HD，ib-sd
nblock =1; % #of fde blocks in each frame, >1 consider IBI; <=1 without considering IBI
Len_chan=30;
M=3; % Alphabet size
Num_chan=10000;
SER_MMSE=zeros(Num_chan,M,length(SNR));
SER_MMSE_mean=zeros(Num_chan,M,length(SNR));
BER_target=zeros(Num_chan*M*length(SNR),1);
Input=zeros(Num_chan,Len_chan,M,length(SNR));
train_data=zeros(Num_chan*M*length(SNR),Len_chan*2+2);
N_SNR=length(SNR);
for index_chan=1:Num_chan
    
    rng(index_chan)
    taps=randi([1,Len_chan],1,1);
    chan=Generat_Channel(taps,Len_chan,index_chan);
    chan_add=zeros(1,chan_order);
    chan_add(1:length(chan))=chan;
    for mod=1:M
        M_mod=2^mod;
        
        for n = 1:length(SNR)
            Input(index_chan,1:Len_chan,mod,n)=chan;
            index=(index_chan-1)*M*N_SNR+(mod-1)*M+n
            train_data(index,1)=n;
            train_data(index,2)=mod;
            train_data(index,3:2+Len_chan)=real(chan);
            train_data(index,2+Len_chan+1:2+Len_chan*2)=imag(chan);
            for m = 1:numTrials % quasi-static channel modeling
                
                
                
                x = randi([0,M_mod-1], N_data, nblock);
                xtx_pload =pskmod(x,M_mod) ;
                %         xtx_pload=ones(128,100);
                xtx_pload_Power = sum(abs(xtx_pload).^2,1)/size(xtx_pload,1); % check transmitted signal power along each tx antenna
                
                PN_add=[];
                for i=1:nblock
                    PN_add=[PN_add CP.'];
                end
                xtx = [xtx_pload;PN_add];  %% add PN
                
                %%
                xtx = [CP.';  reshape(xtx, [], 1)];
                
                fadesig = filter(chan,1, xtx); % Effect of channel, quasi-static
                
                %%
                xtxPower = sum(abs(xtx(:)).^2)/length(xtx(:)); % check transmitted signal power
                xtxPower_dB = 10*log10(xtxPower); % convert to dB
                fadesigPower = sum(abs(fadesig(:)).^2)/length(fadesig(:)); % check fading signal power
                fadesigPower_dB = 10*log10(fadesigPower); % convert to dB
                %         noisePower_dB = fadesigPower_dB-SNR(n); % compute noise power, in dB
                noisePower_dB = xtxPower_dB-SNR(n);       % compute noise power, with
                noisePower = 10^(noisePower_dB/10); % convert backto linear scale
                rng(numTrials)
                noise=sqrt(noisePower)*sqrt(1/2)*(randn(size(fadesig))+1i*randn(size(fadesig)));
                %         rnoise = awgn(fadesig,SNR(n),'measured'); % Additive receiver noise, receiver snr
                %         rnoise = awgn(fadesig, SNR(n), xtxPower_dB); % additive noise with the assumption
                rnoise  =  fadesig+ noise;
                
                
                h_CE=CE_IPNLMS(rnoise(1:N_cp),CP,chan_order);
                
                mse1=sum(abs(h_CE-chan_add).^2)/sum(abs(chan_add).^2);
                %
                %         Gamma=1;lammda=0.1;
                %         h_CE_DR=CE_IPNLMS_data_reuse(rnoise(1:N_cp),CP,chan_order,Gamma,lammda);
                %         mse2=sum(abs(h_CE_DR-chan_add).^2)/sum(abs(chan_add).^2)
                %
                %         Gamma=5;lammda=0.9;
                %         h_CE_DR=CE_IPNLMS_data_reuse(rnoise(1:N_cp),CP,chan_order,Gamma,lammda);
                %         mse3=sum(abs(h_CE_DR-chan_add).^2)/sum(abs(chan_add).^2)
                
                %             h_CE=h_CE;
                rnoise = rnoise(N_cp+1:end, :); % remove PN
                rnoise_CP = reshape(rnoise, [], nblock);  %%
                
                
                
                
%                 H = fft(chan, nfft_PN, 2);
                H=fft(h_CE, nfft_PN, 2);
                
                for i_block=1:nblock
                    
                    rnoise_inter=rnoise_CP(:,i_block);  %% 取一个包
                    fdein = fft(rnoise_inter,nfft_PN); % convert into frequency domain for FDE
                    
                    H_mmse = (H')./((abs(H).^2+ones(size(H))/(xtxPower/noisePower)).');% mmse with transmitted signal to receiver noise power ratio
                    fdeout_mmse = fdein.*H_mmse;
                    xrx_mmse = ifft(fdeout_mmse,nfft_PN); % convert back to time domain for detection
                    %                 z_mmse=demodulate(modem.qamdemod(M_mod),xrx_mmse);
                    z_mmse = pskdemod(xrx_mmse,M_mod);
                    SER_MMSE_temp=BER_Cacula(z_mmse(1:N_data),x(:,i_block));
                    SER_MMSE(index_chan,mod,n) =SER_MMSE(index_chan,mod,n) +SER_MMSE_temp;
                    
                    %                 H_mmse_CE = (H_CE')./((abs(H_CE).^2+ones(size(H_CE))/(xtxPower/noisePower)).');% mmse with transmitted signal to receiver noise power ratio
                    %                 fdeout_mmse_CE = fdein.*H_mmse_CE;
                    %                 xrx_mmse_CE = ifft(fdeout_mmse_CE,nfft_PN); % convert back to time domain for detection
                    %                 z_mmse_CE=demodulate(modem.qamdemod(M_mod),xrx_mmse_CE);
                    %
                    %                 SER_MMSE_temp_CE=BER_Cacula(z_mmse_CE(1:N_data),x(:,i_block));
                    %                 SER_MMSE_CE(n)=SER_MMSE_CE(n)+SER_MMSE_temp_CE;
                    %
                end
            end
            SER_MMSE_mean(index_chan,mod,n) = SER_MMSE(index_chan,mod,n) /numTrials/nblock;
            %             SER_FDO_MMSE_mean = SER_FDO_MMSE/numTrials/nblock;
            BER_target(index)=SER_MMSE_mean(index_chan,mod,n);
        end
    end
    
    % close(h_waitbar)
    
    %     SER_FDO_MMSE_CE_mean= SER_FDO_MMSE_CE/numTrials/nblock;
    %     SER_MMSE_CE_mean= SER_MMSE_CE/numTrials/nblock;
end
%%

% save SER_MMSE_mean SER_MMSE_mean
save('Data_set_channelEstimation', 'train_data', 'BER_target');
