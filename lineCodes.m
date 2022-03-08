clc;
clear variables;
close all;

A=1;

ntow=3;
Rb=1000;Tb=1/Rb;
noOfBits=200;

fsample=1000*Rb;Tsample=1/fsample;

t=0:Tsample:noOfBits*Tb-Tsample;
n=0:noOfBits-1;
delay=0:Tb:Tb*noOfBits-Tb;


alpha = 1;                     %For Raised Cosine Pulse


%bit generation
bitStream=randi([0,1],1,noOfBits);

% %Unipolar NRZ
% %basic pulse and scalar
% basicPulse=rectpuls(t-Tb/2,Tb);
% scalingFactor=A*(bitStream==1);

% %Polar NRZ
% %basic pulse and scalar
% basicPulse=rectpuls(t-Tb/2,Tb);
% scalingFactor=A*(bitStream==1)-A*(bitStream==0);

% %Unipolar RZ
% %basic pulse and scalar
% basicPulse=rectpuls(t-Tb/2,Tb/2);
% scalingFactor=A*(bitStream==1);

% %Polar RZ
% %basic pulse and scalar
% basicPulse=rectpuls(t-Tb/2,Tb/2);
% scalingFactor=A*(bitStream==1)-A*(bitStream==0);

% %Manchester
% %basic pulse and scalar
% basicPulse=rectpuls(t-Tb/4,Tb/2)-rectpuls(t-3*Tb/4,Tb/2);
% scalingFactor=A*(bitStream==1)-A*(bitStream==0);

% % %Ideal Nyquist pulse
% tt=[-ntow*Tb:Tsample:ntow*Tb-Tsample];
% basicPulse=sinc(Rb*tt);   %only using 2 zero crossings of sinc pulse
% basicPulse(numel(t))=0;
% scalingFactor=A*(bitStream==1)-A*(bitStream==0);

%Raised Cosine Pulse
tt=-ntow*Tb : Tsample : ntow*Tb-Tsample;
basicPulse=zeros(1,length(tt));
for i=1:length(tt)
    if(tt(i)==Tb/(2*alpha) || tt(i)== -Tb/(2*alpha)) %when t=Tb/(2*alpha) the denominator term
                                                     %term is zero. This is
                                                     %the limit of the
                                                     %funciton as t
                                                     %approaches
                                                     %Tb/(2*alpha)
        basicPulse(i)= (pi/4)*sinc(1/(2*alpha));
    else
        basicPulse(i)= (sinc(tt(i)/Tb)*cos(pi*alpha*tt(i)/Tb))/(1-(4*alpha^2*tt(i)^2)/(Tb^2));
    end
end

basicPulse(numel(t))=0;                %To make the size of basic pulse large enough
                                       %so that pulstran can work properly
scalingFactor=A*(bitStream==1)-A*(bitStream==0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

train=pulstran(t,[delay;scalingFactor]',basicPulse,fsample);
[ps,f]=pspectrum(train,fsample,'FrequencyResolution',100);

figure;
sgtitle("Time domain representation");
subplot(3,1,1);
plot(t,basicPulse);
title('Basic pulse');
xlabel('time (in s)');ylabel('amplitude');
axis([0,noOfBits*Tb,-1.5,1.5]);

subplot(3,1,2);
stem(n,bitStream);
title('Bit stream');
xlabel('index');ylabel('binary value');

subplot(3,1,3);
%plot(t,train);
plot(t,train,t,t*0+1,t,t*0-1);
title('Line coded signal');
xlabel('time (in s)');ylabel('amplitude');
axis([0,noOfBits*Tb,-1.5,1.5]);

figure;
sgtitle("Power Spectrum for random binary input");
plot(f/1000,10*log10(ps));
xlabel('Frequency(in kHz)');ylabel('Power(in dB)');
axis([0,10,-50,3]);