clc;
clear all;

%Frequency initialization
fc=200e3;
fm=5e3; 

fsample=16*fc;

%Amplitude initialization
Am=0.25;
Ac=sqrt(2);

%Frequency sensitivity, Frequency deviation and Beta
kf=100e3;

delf = kf*Am;
beta = delf/fm;
% beta = 0.1

%time vector
t=0:1/fsample:1000e-3;


ct=Ac*cos(2*pi*fc*t);
mt=Am*cos(2*pi*fm*t);

st=Ac*cos(2*pi*fc*t + beta*sin(2*pi*fm*t));

rt=fmdemod(st,fc,fsample,delf);
display(beta);

figure('Name','Time domain analysis');
subplot(4,1,1);
plot(t,mt);
title("Message Signal");
xlabel("time [in s]");ylabel("m(t)");
xlim([0 400e-6]);

subplot(4,1,2);
plot(t,ct);
title("Carrier Signal");
xlabel("time [in s]");ylabel("c(t)");
xlim([0 400e-6]);

subplot(4,1,3);
plot(t,st);
title("FM Modulated Signal");
xlabel("time [in s]");ylabel("s(t)");
xlim([0 400e-6]);

subplot(4,1,4);
plot(t,rt);
title("FM Demodulated Signal");
xlabel("time [in s]");ylabel("r(t)");
xlim([0 400e-6]);

figure('Name','Power Spectrum');
[ps,f]=pspectrum(st,fsample,'FrequencyResolution',50);
plot(f,10*log10(ps));
title("Power spectrum of FM signal");
xlabel('frequency [in Hz]');
ylabel('Pfm(f) [in dB]');
hold on;
plot(f,(f*0)-40)
grid on;
xlim([100e3 300e3]);
