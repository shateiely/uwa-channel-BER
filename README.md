# uwa-channel-BER
运行BER_data_generate_new.m
第111行：H=fft(chan, nfft_PN, 2); 为生成已知信道
第112行：H = fft（h_CE，nfft_PN，2）; 为生成未知信道
生成表格train_data第一列为调制阶数，其余为信道数据
BER_target为对应调制阶数和信道下的误码率
