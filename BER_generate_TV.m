% �����������ƺ��ŵ��Ĺ�ϵ������
% �ŵ�ʹ������ŵ�
% ���Ʒ�ʽΪ BPSK QPSK 8PSK
% ͨ������Ϊ���ز�Ƶ�����
% ���ն����ȼ����ŵ���֪

clc
close all
clear

Tc=0.25*1e-3; % ���ż��

Len_block = 1024; % ����CP�����ݵ��ܳ���

chan_order=50; % ���Ƶ��ŵ����� 
N_cp=256;

PN=randi([0,1], 1, 512); % 1x512��ȡֵ0-1��α�������-->��512��α�����
temp=[PN(1,1:N_cp)   ]; % ȡPN��ǰһ��(256)
CP=temp+1j*temp; % 1->1+1j,0->0+0j
N_cp =length(CP); % cp�ĳ���
N_data = Len_block-N_cp; % ���ݳ���=CP�����ݵ��ܳ���-CP����

NumTrials = 1; % #of frames in each evaluation

nfft1=N_data;          % ԭʼ���ݿ鳤�� 
nfft_PN=nfft1+N_cp;   % ����PN�����fft����


% 
SNR = 12:2:16; % ���������12��14��16
SER = zeros(4,length(SNR)); % zf, mmse, dfe, ib-HD，ib-sd
nblock =10; % >1����IBI; <=1��IBI
Len_chan=30; % �ŵ��ྶ���� 1:30
M=3; % ���ƽ���
Num_chan=10000;
SER_MMSE=zeros(Num_chan,M,length(SNR));
SER_MMSE_mean=zeros(Num_chan,M,length(SNR));

N_SNR=length(SNR);

rng('default') % ��ʼ���������������ʹ����߱����ظ��ԡ�
N_data_pilot = round((Len_block-N_cp)/M); % ��Ƶ����768/M3=256
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
    taps=randi([1,Len_chan],1,1);
    chan=Generat_Channel(taps,Len_chan,index_chan);   
    Delta_Doppler_taps = rand(1)-0.5;  % ������Ƶ��
    
    for mod=1:M
        M_mod=2^mod; % ���ֵ��Ʒ�ʽ
        for n = 1:length(SNR) % ��������� 14 16 18
            for m = 1:NumTrials % quasi-static channel modeling
                
                x = randi([0,M_mod-1], N_data, nblock);
                xtx_pload =pskmod(x,M_mod) ;
                %         xtx_pload=ones(128,100);
                xtx_pload_Power = sum(abs(xtx_pload).^2,1)/size(xtx_pload,1); % check transmitted signal power along each tx antenna
                
                PN_add=[];
                for i=1:nblock
                    PN_add=[PN_add CP.']; % CP:256x1,nblock��CP
                end
                xtx = [xtx_pload;PN_add];  %% add PN
                xtx_temp= reshape(xtx,[],1);
                xtx_all= [Pilot;xtx_temp];
                index=(index_chan-1)*M*N_SNR*NumTrials+(mod-1)*M*NumTrials+(n-1)*NumTrials+m;
                %             train_data(index,1)=n;
                train_data(index,1)=mod;
                % fadesig = filter(chan,1, xtx_all); % Effect of channel, quasi-static
                %% ʱ��
                Channel_TD=zeros(Len_chan,length(xtx_all));
                [N_path,Num_sampling]=size(Channel_TD);
                T_index_temp = repmat(1:Len_block,1,nblock);
                T_index = [[1:Len_x] T_index_temp];
                for itao=1:N_path
                    Channel_TD(itao,:)=chan(itao)*exp(1j*2*pi *T_index*Delta_Doppler_taps*Tc );
                    %             Channel_TD(delay_taps(itao),:)=chan_coef(itao);
                end
                % Channel_TD_norm=Channel_TD;
                %% ��һ��
                for every_time=1:Num_sampling
                    Channel_TD_norm(:,every_time)=Channel_TD(:,every_time)./sqrt(sum(abs(Channel_TD(:,every_time)).^2)); % ��һ��
                end
                %% ���
                Rxtemp=zeros(N_path,N_path+length(xtx_all)-1);
                for itao=1:N_path
                    Rxtemp(itao,itao:length(xtx_all)+itao-1) =  Channel_TD_norm(itao,:).*xtx_all.';
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
                %rnoise  =  fadesig+ noise;
                %rnoise  =  fadesig;
                rnoise  =  fadesig+ noise;
                rnoise_input = rnoise(1:Len_x);
                
                h_CE=CE_IPNLMS(rnoise(1:N_cp),CP,chan_order);
                train_data(index,2:1+length(rnoise_input))=real(rnoise_input); % 估计的信�?
                train_data(index,1+length(rnoise_input)+1:1+length(rnoise_input)*2)=imag(rnoise_input);
                
                
                rnoise = rnoise(Len_x+1:end, :); % remove PN
                rnoise_CP = reshape(rnoise, [], nblock);  %%
                %H=fft(chan, nfft_PN, 2); % ��֪�ŵ�
                H=fft(h_CE, nfft_PN, 2); % δ֪�ŵ�
                
                for i_block=1:nblock
                    
                    rnoise_inter=rnoise_CP(:,i_block);  %% 取一个包
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
save('Data_set_all', 'train_data', 'BER_target');
