clc;
clear variables;

%Frequency initialization
fc=100e3;
fm=5e3;

fsample=16*fc;

%Amplitude initialization
Am=7;
Ac=5;

%Amplitude sensitivity
ka=0.25;

%time vector
t=0:1/fsample:100e-3;


ct=Ac*cos(2*pi*fc*t);
mt=Am*cos(2*pi*fm*t);

st=ct+ka*mt.*ct;

[rt,blabla]=envelope(st);

figure('Name','Time domain analysis');
subplot(4,1,1);
plot(t,mt);
title("Message Signal");
xlabel("time (s)");ylabel("m(t) (V)");

subplot(4,1,2);
plot(t,ct);
title("Carrier Signal");
xlabel("time (s)");ylabel("c(t) (V)");

subplot(4,1,3);
plot(t,st);
title("AM Modulated Signal");
xlabel("time (s)");ylabel("s(t) (V)");

subplot(4,1,4);
plot(t,rt);
title("AM Demodulated Signal");
xlabel("time (s)");ylabel("r(t) (V)");

figure('Name','Power Spectrum');
[ps,f]=pspectrum(st,fsample,'FrequencyResolution',100);
plot(f,10*log10(ps));
title("Power Spectrum of AM signal");
xlabel("frequency(Hz)");ylabel("Pam(f) (dB)");
grid on;
