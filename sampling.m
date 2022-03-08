clc;
clear variables;
close all;

%message signal parameters
fm=2e3; Tm=1/fm;
Am=2;
Adc=2;

%Pulse train parameters
fsample=2*fm; Tsample=1/fsample;
dutyCycle=1/100;
pulseWidth=dutyCycle/fsample;

zeroCrossingFrequency=1/pulseWidth

fsdas=100*fsample;               %Sampling rate used by MATLAB to work with

%time vector
t=0:1/fsdas:500*Tm;


%Message signal generation
x=Adc+Am*cos(2*pi*fm*t);

%Pulse train generation
delays=0:1/fsample:1-1/fsample;
pulseTrain=pulstran(t-pulseWidth/2,delays,'rectpuls',pulseWidth);

%The process of Sampling
v=pulseTrain.*x;


%Filter Design and Filtering of Signal
filterOrder=4;
filterCutOffFreq=fm;
[b,a]= butter(filterOrder,filterCutOffFreq/(fsdas/2));
[H,f]= freqz(b,a,fsdas,fsdas);
y=filter(b,a,v);

%Time domain Results plotting
figure;
sgtitle("Time Domain Analysis of Sampling");
subplot(4,1,1); plot(t,x);
xlabel("time [in s]");ylabel("x(t) [in volt]");
title("Input Analog Signal");
axis([0.1 0.105 min(x) max(x)])

subplot(4,1,2); plot(t,pulseTrain);
xlabel("time [in s]");ylabel("p(t) [in volt]");
title("Pulse train")
axis([0.1 0.105 min(pulseTrain) max(pulseTrain)]);

subplot(4,1,3); plot(t,v);
xlabel("time [in s]");ylabel("s(t) [in volt]");
title("Sampled Signal");
axis([0.1 0.105 min(v) max(v)]);

subplot(4,1,4); plot(t,y);
xlabel("time [in s]");ylabel("r(t) [in volt]");
title("Reconstructed signal");
axis([0.1 0.105 min(y) max(y)])

%Power spectrum estimation of the signals
fres=fm/20;
[px,fx]= pspectrum(x,fsdas,'FrequencyResolution',fres);
[pv,fv]= pspectrum(v,fsdas,'FrequencyResolution',fres);
[py,fy]= pspectrum(y,fsdas,'FrequencyResolution',fres);

%Frequency domain results plotting
figure;
sgtitle("Frequency Domain analysis of Sampling");
subplot(4,1,1); plot(fx,10*log10(px));
xlabel("frequency [in Hz]");ylabel("X(f)");
title("Input signal's spectum");
axis([0 14e4 -100 10]);

subplot(4,1,3); plot(fv,10*log10(pv));
xlabel("frequency [in Hz]");ylabel("S(f)");
title("Sampled signal's spectrum");
axis([0 14e4 -100 10]);

subplot(4,1,2); semilogx(f,10*log10(abs(H)));
xlabel("frequency [in Hz]");ylabel("H(f)");
title("Filter's Frequency Response");
axis([0 14e4 -100 10]);

subplot(4,1,4); plot(fy,10*log10(py));
xlabel("frequency [in Hz]");ylabel("R(f)");
title("Reconstructed signal's spectrum");
axis([0 14e4 -100 10]);