function Channel_Real=Choose_channel(Len_chan,index_chan,xtx_all)
if index_chan <= 4
    load(strcat('C:\Users\dcfor\Documents\MATLAB\uwa-channel&BER\channel\BCH1_',num2str(index_chan),'.mat'), 'h')
    resample_multiple=4; % 抽样倍数
    resample_path=1;
    for h_path=1:size(h,2)
        if rem(h_path,resample_multiple)==0
            h_resample(:,resample_path)=h(:,h_path);
            resample_path=resample_path+1;
        end
    end
    Channel_temp=h_resample(:,1:Len_chan).';
    sampling_temp=size(Channel_temp,2);
    interp_multiple=ceil(length(xtx_all)/sampling_temp); % 插值倍数
    for nn=1:Len_chan
        channel_interp(nn,:) = interp(Channel_temp(nn,:),interp_multiple); %对每列进行multiple倍插值
    end
    Channel_Real_temp=channel_interp(:,1:length(xtx_all)); % 截取与发射信号等长的信道
elseif index_chan >= 5 && index_chan <= 20
    load(strcat('C:\Users\dcfor\Documents\MATLAB\uwa-channel&BER\channel\KAU1_',num2str(index_chan-4),'.mat'), 'h')
    
    resample_multiple=4; % 抽样倍数
    resample_path=1;
    for h_path=1:size(h,2)
        if rem(h_path,resample_multiple)==0
            h_resample(:,resample_path)=h(:,h_path);
            resample_path=resample_path+1;
        end
    end
    Channel_temp=h_resample(:,1:Len_chan).';
    sampling_temp=size(Channel_temp,2);
    interp_multiple=ceil(length(xtx_all)/sampling_temp); % 插值倍数
    for nn=1:Len_chan
        channel_interp(nn,:) = interp(Channel_temp(nn,:),interp_multiple); %对每列进行multiple倍插值
    end
    Channel_Real_temp=channel_interp(:,1:length(xtx_all)); % 截取与发射信号等长的信道
elseif index_chan >= 21 && index_chan <= 36
    load(strcat('C:\Users\dcfor\Documents\MATLAB\uwa-channel&BER\channel\KAU2_',num2str(index_chan-20),'.mat'), 'h')
    
    resample_multiple=4; % 抽样倍数
    resample_path=1;
    for h_path=1:size(h,2)
        if rem(h_path,resample_multiple)==0
            h_resample(:,resample_path)=h(:,h_path);
            resample_path=resample_path+1;
        end
    end
    Channel_temp=h_resample(:,1:Len_chan).';
    sampling_temp=size(Channel_temp,2);
    interp_multiple=ceil(length(xtx_all)/sampling_temp); % 插值倍数
    for nn=1:Len_chan
        channel_interp(nn,:) = interp(Channel_temp(nn,:),interp_multiple); %对每列进行multiple倍插值
    end
    Channel_Real_temp=channel_interp(:,1:length(xtx_all)); % 截取与发射信号等长的信道
elseif index_chan >= 37 && index_chan <= 96
    load(strcat('C:\Users\dcfor\Documents\MATLAB\uwa-channel&BER\channel\NCS1_ (',num2str(index_chan-36),').mat'), 'h')
    
    resample_multiple=4; % 抽样倍数
    resample_path=1;
    for h_path=1:size(h,2)
        if rem(h_path,resample_multiple)==0
            h_resample(:,resample_path)=h(:,h_path);
            resample_path=resample_path+1;
        end
    end
    Channel_temp=h_resample(:,1:Len_chan).';
    sampling_temp=size(Channel_temp,2);
    interp_multiple=ceil(length(xtx_all)/sampling_temp); % 插值倍数
    for nn=1:Len_chan
        channel_interp(nn,:) = interp(Channel_temp(nn,:),interp_multiple); %对每列进行multiple倍插值
    end
    Channel_Real_temp=channel_interp(:,1:length(xtx_all)); % 截取与发射信号等长的信道
elseif index_chan >= 97
    load(strcat('C:\Users\dcfor\Documents\MATLAB\uwa-channel&BER\channel\NOF1_ (',num2str(index_chan-96),').mat'), 'h')
    
    resample_multiple=4; % 抽样倍数
    resample_path=1;
    for h_path=1:size(h,2)
        if rem(h_path,resample_multiple)==0
            h_resample(:,resample_path)=h(:,h_path);
            resample_path=resample_path+1;
        end
    end
    Channel_temp=h_resample(:,1:Len_chan).';
    sampling_temp=size(Channel_temp,2);
    interp_multiple=ceil(length(xtx_all)/sampling_temp); % 插值倍数
    for nn=1:Len_chan
        channel_interp(nn,:) = interp(Channel_temp(nn,:),interp_multiple); %对每列进行multiple倍插值
    end
    Channel_Real_temp=channel_interp(:,1:length(xtx_all)); % 截取与发射信号等长的信道
end

for ii=1:length(xtx_all)
    Channel_Real(:,ii)=Channel_Real_temp(:,ii)./sqrt(sum(abs(Channel_Real_temp(:,ii)).^2)); % 归一化
end
