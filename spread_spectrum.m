% clc;
% clear variables;
% close all;
% 
% Rb = 1e3;               %bit rate
% Nb = 10000;             %no of bits to be generated
% b =2*randi([0 1],Nb,1)-1;   %random binary bits
% b(b==0) = -1;
% 
% PNseq = [0 0 1 1 1 0 1];
% LPN= length(PNseq);
% 
% Rc = LPN*Rb;
% OSRc = 4;
% OSRb = LPN*OSRc;
% 
% fs = OSRb*Rb;
% pulsec = ones(1,OSRc);
% pulseb = ones(1,OSRb);
% PNseq(PNseq==0) = -1;
% chipseq = zeros(1,(LPN - 1)*OSRc +1);
% chipseq(1:OSRc:end) = PNseq;
% chipseq = conv(chipseq,pulsec);
% tc = (0:LPN*OSRc - 1)/(Rc*OSRc);
% 
% bitseq1 = zeros(1,(Nb-1)*OSRb+1);
% bitseq1(1:OSRb:end) = b;
% bitseq_raw = conv(bitseq1,pulseb);
% tb = (0:Nb*OSRb-1)/(Rb*OSRb);
% 
% c=[];
% for ii = 1:Nb
%     c = [c chipseq];
% end
% 
% 
% 
% figure;
% sgtitle('Time Domain Characteristic of DSSS');
% subplot(4,1,1);
% plot(tb,bitseq_raw);
% title('input bitstream');
% xlabel('time (in s)');
% ylabel('amplitude');
% axis([0 max(tb) -1.1 1.1]);
% 
% subplot(4,1,2);
% plot(tb,c); title('pseudorandom noise(PN) signal');
% xlabel('time (in s)');
% ylabel('amplitude');
% axis([0 max(tb) -1.1 1.1]);
% 
% xseq= bitseq_raw.*c;
% y=xseq.*c;
% 
% subplot(4,1,3);
% plot(tb,xseq); title('DSSS - BPSK Signal');
% xlabel('time (in s)');
% ylabel('amplitude');
% axis([0 max(tb) -1.1 1.1]);
% 
% subplot(4,1,4);
% plot(tb,y); title('received signal');
% xlabel('time (in s)');
% ylabel('amplitude');
% axis([0 max(tb) -1.1 1.1]);
% 
% [pb,f] = pwelch(bitseq_raw,[],[],[],fs);
% [px,f] = pwelch(xseq,[],[],[],fs);
% 
% figure;
% sgtitle('Frequency Domain Characteristic of DSSS');
% subplot(2,1,1);
% plot(f,10*log10(pb));
% title('power spectrum of input bit stream');
% xlabel('frequency (in Hz)');
% ylabel('power (in dB)');
% subplot(2,1,2);
% plot(f,10*log10(px));
% title('power spectrumof the spreaded signal');
% xlabel('frequency (in Hz)');
% ylabel('power (in dB)');


clc;
clear variables;
close all;

PN_seq = [0 0 1 1 1 0 1];
LPN = length(PN_seq);

Nb = 1e3;

Rb = 1e3; Tb=1/Rb;
Rc = LPN*Rb; Tc=1/Rc;

OSRc = 4;
Tsample = Tc/OSRc;Fsample=1/Tsample;

t=0:Tsample:Tb*Nb-Tsample;

bk = randi([0 1],1,Nb);
bk_lineCoded = bk;
bk_lineCoded(bk_lineCoded==0) = -1;
bitBasicPulse = rectpuls(t-Tb/2,Tb);
delay = 0:Tb:Nb*Tb-Tb;
bit_genout = pulstran(t,[delay;bk_lineCoded]',bitBasicPulse,Fsample);
% stem(t,bit_genout);

PN_seq_lineCoded = PN_seq;
PN_seq_lineCoded(PN_seq_lineCoded ==0) = -1;

chipBasicPulse = rectpuls(t-Tc/2,Tc);
delay = 0:Tc:1*Tb-Tc;
PNbasicPulse=pulstran(t,[delay;PN_seq_lineCoded]',chipBasicPulse,Fsample);
% stem(t,PNBasicPulse)

delay = 0:Tb:Nb*Tb-Tb;
PN_genout = pulstran(t,delay,PNbasicPulse,Fsample);
% stem(t,PN_genout)

spreaded_seq = PN_genout.*bit_genout;
despreaded_seq = spreaded_seq .* PN_genout;


figure;
sgtitle('Time Domain Characteristic of DSSS');
subplot(4,1,1);
plot(t,bit_genout);
title('input bitstream');
xlabel('time (in s)');
ylabel('amplitude');
axis([0 max(t) -1.1 1.1]);

subplot(4,1,2);
plot(t,PN_genout); title('pseudorandom noise(PN) signal');
xlabel('time (in s)');
ylabel('amplitude');
axis([0 max(t) -1.1 1.1]);

subplot(4,1,3);
plot(t,spreaded_seq); title('DSSS - BPSK Signal');
xlabel('time (in s)');
ylabel('amplitude');
axis([0 max(t) -1.1 1.1]);

subplot(4,1,4);
plot(t,despreaded_seq); title('received signal');
xlabel('time (in s)');
ylabel('amplitude');
axis([0 max(t) -1.1 1.1]);

% [pb,f] = pwelch(bit_genout,[],[],[],Fsample);
% [px,f] = pwelch(spreaded_seq,[],[],[],Fsample);

[pb,f] = pspectrum(bit_genout,Fsample,'FrequencyResolution',5);
[px,f] = pspectrum(spreaded_seq,Fsample,'FrequencyResolution',5);


figure;
sgtitle('Frequency Domain Characteristic of DSSS');
subplot(2,1,1);
plot(f,10*log10(pb));
title('power spectrum of input bit stream');
xlabel('frequency (in Hz)');
ylabel('power (in dB)');
subplot(2,1,2);
plot(f,10*log10(px));
title('power spectrumof the spreaded signal');
xlabel('frequency (in Hz)');
ylabel('power (in dB)');
