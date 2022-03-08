clc;
clear variables;
close all;

M=2;               %number of symbols
N=1e6;             %N=1e3 for signal constellation
n = log2(M);       %no. of bits / symbol

bk = randi([0 1],N*n,1);   %bk - bit stream

dpskmod = comm.DPSKModulator(2,pi,'BitInput',true);
dpskdemod = comm.DPSKDemodulator(2,pi,'BitOutput',true);

xs = dpskmod(bk);             %xs - electrical value
scatter(real(xs),imag(xs),'bo','filled');
title('Constellation diagram of DPSK');
xlabel('in-phase'); ylabel('quadrature-phase');

%adding error
SNRdB=1:1:15;
Eb_by_No_dB=zeros(1,length(SNRdB));
BER = zeros(1,length(SNRdB));
BER_th = zeros(1,length(SNRdB));

for i=1:length(SNRdB)
    SNRi = SNRdB(i);
    rt = awgn(xs,SNRi,'measured');
    
%     figure;
%     plot(real(rt),imag(rt),'b.',real(xs),imag(xs),'r*');  
%                                                     %plotting where the corrupted
%                                                     %symbol lands up
%     legend('received','transmitted');
%     title("DPSK constellation at SNR="+num2str(SNRi)+" dB");
%     axis([-3.5, 3.5, -1.5, 1.5]);
%     xlabel('real');ylabel('imaginary');
    
    %br stands for bit recived
    br=dpskdemod(rt);  
    snr = 10^(SNRi/10);
    Eb_by_No = snr;
    Eb_by_No_dB(i) = 10*log10(Eb_by_No);
    
    BER_th(i) = 0.5*exp(-Eb_by_No);
    BER(i) = length(find(bk ~= br))/N;   
end

figure;
semilogy(Eb_by_No_dB,BER,'b-',Eb_by_No_dB,BER_th,'r*');
legend('BER from simulation','BER from theory');
axis([0, 12, 10^-7 1]);grid on;
xlabel('Eb/No (in dB)'); ylabel('bit error probability');
title('BER performance of DPSK under AWGN');


