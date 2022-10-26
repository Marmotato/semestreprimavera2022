
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script for plotting the time domain and frequency domain representation
% of raised cosine filters for various values of alpha
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
fs =10;
NFFT=1024;
FPulse=[-NFFT/2:NFFT/2-1]/NFFT;
alpha = 0.5;
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%defining the RC filter
sincNum = sin(pi*[-fs:1/fs:fs]); % numerator of the sinc function
sincDen = (pi*[-fs:1/fs:fs]); % denominator of the sinc function
sincDenZero = find(abs(sincDen) < 10^-10);
sincOp = sincNum./sincDen;
sincOp(sincDenZero) = 1; % sin(pix/(pix) =1 for x =0

cosNum = cos(alpha*pi*[-fs:1/fs:fs]);
cosDen = (1-(2*alpha*[-fs:1/fs:fs]).^2);
cosDenZero = find(abs(cosDen)<10^-10);
cosOp = cosNum./cosDen;
cosOp(cosDenZero) = pi/4;

RC_Filter = sincOp.*cosOp;
RC_Filter_Spectrum = (abs((fft(RC_Filter,NFFT))/fs));
RC_Filter_Spectrum_dB = 20*log10((fft(RC_Filter,1024))/fs);
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%LINEAR COMBINATION PULSE
%alpha*PLP(n=1) +(1-alpha)RC
constante=1; %alpha

Exp1_Num=((1-constante)*sin(pi*[-fs:1/fs:fs]).*cos(pi*alpha*[-fs:1/fs:fs])); %RC
Exp1_Den=pi*[-fs:1/fs:fs].*(1-(2*alpha*[-fs:1/fs:fs]).^2);
Exp1DenZero=find(abs(Exp1_Den) < 10^-10);
Exp1_Total=Exp1_Num./Exp1_Den;
Exp1_Total(Exp1DenZero)=0.5;

Exp2_Num=constante*sin(pi*[-fs:1/fs:fs]).*sin(pi*alpha*[-fs:1/fs:fs]); %PLP
Exp2_Den=pi^2*alpha*[-fs:1/fs:fs].^2;
Exp2DenZero=find(abs(Exp2_Den) < 10^-10);
Exp2_Total=Exp2_Num./Exp2_Den;
Exp2_Total(Exp2DenZero)=0.5;

Linear_Combination_Pulse_PLP=(Exp1_Total+Exp2_Total);
Linear_Combination_Pulse_PLP_Spectrum=(abs((fft(Linear_Combination_Pulse_PLP,NFFT))/fs));
Linear_Combination_Pulse_PLP_Spectrum_dB = 20*log10((fft(Linear_Combination_Pulse_PLP,1024))/fs);



%EXPONENTIAL LINEAR PULSE
constante=1; %alpha
beta = 0.1 %Beta


Exp2_Num=constante*sin(pi*[-fs:1/fs:fs]).*sin(pi*alpha*[-fs:1/fs:fs]); %PLP
Exp2_Den=pi^2*alpha*[-fs:1/fs:fs].^2;
Exp2DenZero=find(abs(Exp2_Den) < 10^-10);
Exp2_Total=Exp2_Num./Exp2_Den;
Exp2_Total(Exp2DenZero)=0.5;

Exp3 = exp(-pi*(beta/2)*[-fs:1/fs:fs])



Exponential_Linear_Pulse=( Exp3.*Exp2_Total);
Exponential_Linear_Pulse_Spectrum=(abs((fft(Linear_Combination_Pulse_PLP,NFFT))/fs));
Exponential_Linear_Pulse_Spectrum_dB = 20*log10((fft(Linear_Combination_Pulse_PLP,1024))/fs);




%BETTER THAN RAISED COSINE PULSE

beta2 = (pi * alpha) / log(2);

% Calculate the sinc expression
sinc_BTRC = sin([-fs:1/fs:fs])./[-fs:1/fs:fs];
sinc_BTRC_zero=find(abs(sinc_BTRC) < 10^-10);

% Calculate the second expression for BTRC
numerador_second_BTRC = (2 * beta2 * [-fs:1/fs:fs]) .* sin(pi * alpha * [-fs:1/fs:fs]);
numerador_second_BTRC = numerador_second_BTRC + 2 .* cos(pi * alpha * [-fs:1/fs:fs]) - 1;

% Never 0
denominador_second_BTRC = (beta2 * [-fs:1/fs:fs]).^2 + 1;

% Get the BTRC Pulse
BTRC = sinc_BTRC .* numerador_second_BTRC ./ denominador_second_BTRC;
BTRC(sinc_BTRC_zero) = numerador_second_BTRC / denominador_second_BTRC;
BTRC_Pulse_Spectrum=(abs((fft(BTRC,NFFT))/fs));
BTRC_Pulse_Spectrum_dB = 20*log10((fft(BTRC,1024))/fs);



%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Respuesta al Impulso
close all
figure(1)
plot([-fs:1/fs:fs],[RC_Filter],'--kd','LineWidth',1)
hold on
plot([-fs:1/fs:fs],[Linear_Combination_Pulse_PLP],'--ko','LineWidth',1)
hold on
plot([-fs:1/fs:fs],[BTRC],'--ks','LineWidth',1)
hold on
plot([-fs:1/fs:fs],[Exponential_Linear_Pulse],'--kp','LineWidth',1)
hold on

title('Impulse Response for ISI-Free Pulses with alpha='+string(alpha),'LineWidth',.1);
legend('RC', 'Linear Combination Pulse', 'BTRC', 'Exponential Linear Pulse PLP');
axis([0 3 -.2 1.1])
hold all

grid on
xlabel('Time, t/T')
ylabel('Amplitude, h(t)')
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Respuesta en Frecuencia
figure(2)
plot(FPulse*2*fs, fftshift(RC_Filter_Spectrum),'--kd','LineWidth',1);
hold on
plot(FPulse*fs*2, fftshift(Linear_Combination_Pulse_PLP_Spectrum),'--ko','LineWidth',1);
hold on
plot(FPulse*fs*2, fftshift(BTRC_Pulse_Spectrum),'--ks','LineWidth',1);
hold on
plot(FPulse*fs*2, fftshift(Exponential_Linear_Pulse_Spectrum),'--kp','LineWidth',1);
hold on
title('Frequency Response for ISI-Free Pulses with alpha='+string(alpha),'LineWidth',.1);
legend('RC', 'Linear Combination Pulse', 'BTRC', 'Exponential Linear Pulse PLP');
axis([-1.6 1.6 0 1])

grid on
xlabel('f/B')
ylabel('H(f)/T')
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% %Spectral Magnitude
% figure(3)
% plot(FPulse*2*fs, fftshift(RC_Filter_Spectrum_dB),'--kd','LineWidth',1);
% hold on
% plot(FPulse*fs*2, fftshift(Linear_Combination_Pulse_PLP_Spectrum_dB),'--ko','LineWidth',1);
% hold on
% 
% legend('RC','Linear Combination Pulse');
% 
% 
% grid on
% xlabel('f/B')
% ylabel('dB')
%--------------------------------------------------------------------------
%-------------------------------------------------------------------------
% Welch spectral estimator.
%RAISED COSINE
% figure(4)
% Fs=10;
% h = spectrum.welch('Hamming') ;    %Create a Welch spectral estimator using a Hamming Window
% Hpsd_RC = psd(h,RC_Filter,'Fs',Fs','CenterDC',true);   %Calculate the PSD 
% normalizefreq(Hpsd_RC); 
% plot (Hpsd_RC)
% hold on
% grid on
% legend('Raised Cosine Pulse');
% 
% figure(5)
% %Linear Combination Pulse
% h = spectrum.welch('Hamming');   %Create a Welch spectral estimator using a Hamming Window
% Hpsd_Linear_Combination_Pulse_PLP = psd(h,Linear_Combination_Pulse_PLP,'Fs',Fs','CenterDC',true);  %Calculate the PSD 
% normalizefreq(Hpsd_Linear_Combination_Pulse_PLP); 
% plot(Hpsd_Linear_Combination_Pulse_PLP)
% hold on
% grid on
% legend('Linear Combination Pulse');


function [y]=sinc_m(x)
    if x==0
        y=[1]
    else
        y = (sin([x])/[x])
        
    end
end