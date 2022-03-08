clc;
close all;

fm = 1e3;               %input message signal frequency
fsample = 1024*fm;      %Sampling rate of signal
Am = 1;                 %message signal amplitude

n=[15];                %Testing 1, 2 and 3 bits

t=0:1/fsample:2/fm-1/fsample;   %To observe 2 cycles of input signal

x=Am*cos(2*pi*fm*t);    %input messag signal
x_power=x*x';

SQNRdB=zeros(1,length(n));      %To store practical SQNR
SQNRdB_th=zeros(1,length(n));   %To store theoretical SQNR

for i=1:length(n)
    delta = 2*Am/(2^n(i)-1);    %Step size of the quantizer
    
    %Mid-tread quantizer
    partition=-Am+delta/2:delta:Am-delta/2;
    codebook=-Am:delta:Am;
    
%     %Mid-rise quantizer
%     partition=-Am+delta:delta:Am-delta;
%     codebook=-Am+delta/2:delta:Am-delta/2;
    
    [index,xq] = quantiz(x,partition,codebook);
    qe=x-xq;
    qe_avg=mean(qe);
    qe_power=qe*qe';
    
    SQNR = x_power/qe_power;
    SQNRdB(i)=10*log10(SQNR);
    SQNRdB_th(i)=6*n(i)+1.72;
    
    subplot(length(n),1,i);
    plot(t,x,'r--',t,xq,'b-');      %plot input signal and
                                    %its quantized form
    title("Time Domain Signal for n="+num2str(n(i)));
    xlabel('time(s)'); ylabel('amplitude');
    legend('message signal','quantized signal');
end

% figure;
% plot(n,SQNRdB,'r-*',n,SQNRdB_th,'b-*');
% title('SQNR');
% xlabel('No of bits used'); ylabel('SQNR(dB)');
% legend('SQNR','SQNR theoretical');


