clc;
clear variables;
close all;

M=4;               %number of symbols
N=1e6;             %N=1e3 for signal constellation
n = log2(M);       %no. of bits / symbol

bk = randi([0 1], 1,N*n);   %bk - bit stream
xs = bk; xs(xs==0)=-1;      %xs - electrical value
                            %Performing polar NRZ

xs_odd = xs(1:2:end);
xs_even = xs(2:2:end);
xs_complex = xs_odd + 1j*xs_even;

scatter(real(xs_complex),imag(xs_complex),'bo','filled');
title('Constellation diagram of QPSK');
xlabel('in-phase'); ylabel('quadrature-phase');

%adding error
SNRdB=1:1:20;
Eb_by_No_dB=zeros(1,length(SNRdB));
BER = zeros(1,length(SNRdB));
BER_th = zeros(1,length(SNRdB));

for i=1:length(SNRdB)
    SNRi = SNRdB(i);
    rt = awgn(xs_complex,SNRi,'measured');
    
%     figure;
%     plot(real(rt),imag(rt),'b.',real(xs_complex),imag(xs_complex),'r*');  
%                                                     %plotting where the corrupted
%                                                     %symbol lands up
%     legend('received','transmitted');
%     title("QPSK constellation at SNR="+num2str(SNRi)+" dB");
%     axis([-3.5, 3.5, -1.5, 1.5]);
%     xlabel('real');ylabel('imaginary');
    
    %br stands for bit recived
    br_odd = real(rt); br_odd(br_odd>=0)=1; br_odd(br_odd<0)=0;         %ML decoding
    br_even = imag(rt); br_even(br_even>=0)=1;br_even(br_even<0)=0;
    
    br=zeros(1,length(xs));
    br(1:2:end)=br_odd;
    br(2:2:end)=br_even;
    
    snr = 10^(SNRi/10);
    Eb_by_No = snr/2;
    Eb_by_No_dB(i) = 10*log10(Eb_by_No);
    
    BER_th(i) = 1*erfc(sqrt(Eb_by_No));
   %BER_th(i) = 0.5*(erfc(sqrt(Eb_by_No))-0.25*(erfc(sqrt(Eb_by_No))^2));
    BER(i) = length(find(bk ~= br))/N;   
end

figure;
semilogy(Eb_by_No_dB,BER,'b-',Eb_by_No_dB,BER_th,'r*');
legend('BER from simulation','BER from theory');
axis([0, 12, 10^-7 1]);grid on;
xlabel('Eb/No (in dB)'); ylabel('bit error probability');
title('BER performance of QPSK under AWGN');