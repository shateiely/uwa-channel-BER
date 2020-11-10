clc
clear
close all
load Data_Predict.mat
[P,q]=size(Real_BER);
Nonzero = find( Real_BER);
Zero_index=setdiff( [1:P],Nonzero) ;
% Real_BER(Zero_index) = 1e-5;
x=1:5000;
semilogy(Real_BER(x),'-bo','LineWidth',1,'MarkerSize',4)
hold on
semilogy(Predict_BER(x),'-r^','LineWidth',1,'MarkerSize',4)
xlabel('次数')
ylabel('误码率')
legend('Real BER','Predict BER')
grid on
figure
Bais=(Predict_BER-Real_BER);
plot(Bais)
MAPE = mean(abs(Bais)./Real_BER)/100