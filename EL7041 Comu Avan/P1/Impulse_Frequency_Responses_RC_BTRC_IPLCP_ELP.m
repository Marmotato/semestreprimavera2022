
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script for plotting the time domain and frequency domain representation
% of raised cosine filters for various values of alpha
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
fs =10;
NFFT=1024;
FPulse=[-NFFT/2:NFFT/2-1]/NFFT;
alpha = 0.5; %roll-off factor, reemplazar para evaluar

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%% RAISED COSINE

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

%% IMPROVED LINEAR COMBINATION PULSE
% ((mu*A +(1-mu)B) * sinc * exp)^gamma

mu = 1.60;
eps = 0.1;
gamma = 1;

sincNum = sin(pi*[-fs:1/fs:fs]); % numerator of the sinc function
sincDen = (pi*[-fs:1/fs:fs]); % denominator of the sinc function
sincDenZero = find(abs(sincDen) < 10^-10);
sincOp = sincNum./sincDen;
sincOp(sincDenZero) = 1; % sin(pix/(pix) =1 for x =0

ExpIPLC = exp(-eps*pi^2*([-fs:1/fs:fs]).^2);

IPLCNum = 4*(1 - mu)*(sin(alpha*pi*[-fs:1/fs:fs]/2).^2) + alpha*pi*mu*[-fs:1/fs:fs].*sin(alpha*pi*[-fs:1/fs:fs]);
IPLCDen = (alpha*pi*[-fs:1/fs:fs]).^2;
IPLCDenZeros = find(abs(IPLCDen) < 10^-10);
IPLC_1 = IPLCNum ./ IPLCDen;
IPLC_1(IPLCDenZeros) = 1;

IPLC = ExpIPLC.*(sincOp.*IPLC_1).^gamma;
IPLC_Spectrum=(abs((fft(IPLC,NFFT))/fs));
IPLC_Spectrum_dB = 20*log10((fft(IPLC,1024))/fs);

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%% EXPONENTIAL LINEAR PULSE
constante=1; %alpha
beta = 0.1; %Beta

Exp2_2_Num=constante*sin(pi*[-fs:1/fs:fs]).*sin(pi*alpha*[-fs:1/fs:fs]); %PLP
Exp2_2_Den=pi^2*alpha*[-fs:1/fs:fs].^2;
Exp2_2DenZero=find(abs(Exp2_2_Den) < 10^-10);
Exp2_2_Total=Exp2_2_Num./Exp2_2_Den;
Exp2_Total(Exp2_2DenZero)=1;

Exp3 = exp(-pi*(beta/2)*[-fs:1/fs:fs].^2);

Exp3_Zero=find(abs(Exp3) < 10^-10);
Exp3(Exp3_Zero) = 1;

Exponential_Linear_Pulse=( Exp3.*Exp2_2_Total);
Exponential_Linear_Pulse(isnan(Exponential_Linear_Pulse))=1;
Exponential_Linear_Pulse_Spectrum=(abs((fft(Exponential_Linear_Pulse,NFFT))/fs));
Exponential_Linear_Pulse_Spectrum_dB = 20*log10((fft(Exponential_Linear_Pulse,1024))/fs);

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%% BETTER THAN RAISED COSINE PULSE
% 
beta2 = (pi * alpha) / log(2);

% Calculate the sinc expression
sinc_BTRC_num = sin(pi.*[-fs:1/fs:fs]);
sinc_BTRC_den = pi.*[-fs:1/fs:fs];
sinc_BTRC_tot = sinc_BTRC_num./sinc_BTRC_den;
sinc_BTRC_den_zero=find(abs(sinc_BTRC_tot) < 10^-10);
sinc_BTRC_tot(sinc_BTRC_den_zero) = 0;

% Calculate the second expression for BTRC
numerador_second_BTRC = (2 * beta2 * [-fs:1/fs:fs]) .* sin(pi * alpha * [-fs:1/fs:fs]);
numerador_second_BTRC = numerador_second_BTRC + (2 .* cos(pi * alpha * [-fs:1/fs:fs]) - 1);

% Never 0
denominador_second_BTRC = (beta2 * [-fs:1/fs:fs]).^2 + 1;

% Get the BTRC Pulse
BTRC = ((sinc_BTRC_tot .* numerador_second_BTRC) ./ denominador_second_BTRC);
%BTRC = arrayfun(@(x) BTRC(x,alpha),[-fs:1/fs:fs]);

BTRC(isnan(BTRC))=1; %corregimos en el punto 0

BTRC_Pulse_Spectrum=(abs((fft(BTRC,NFFT))/fs));
BTRC_Pulse_Spectrum_dB = 20*log10((fft(BTRC,1024))/fs);



%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%% Respuesta al Impulso
close all
figure(1)
plot([-fs:1/fs:fs],[RC_Filter],'--kd','LineWidth',1)
hold on
plot([-fs:1/fs:fs],[IPLC],'--ko','LineWidth',1)
hold on
plot([-fs:1/fs:fs],[BTRC],'--ks','LineWidth',1)
hold on
plot([-fs:1/fs:fs],[Exponential_Linear_Pulse],'--kp','LineWidth',1)
hold on

title('Impulse Response for ISI-Free Pulses with alpha='+string(alpha),'LineWidth',.1);
legend('RC', 'IPLCP', 'BTRC', 'ELP');
axis([0 3 -.2 1.1])
hold all

grid on
xlabel('Time, t/T')
ylabel('Amplitude, h(t)')
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%% Respuesta en Frecuencia
figure(2)
plot(FPulse*2*fs, fftshift(RC_Filter_Spectrum),'--kd','LineWidth',1);
hold on
plot(FPulse*fs*2, fftshift(IPLC_Spectrum),'--ko','LineWidth',1);
hold on
plot(FPulse*fs*2, fftshift(BTRC_Pulse_Spectrum),'--ks','LineWidth',1);
hold on
plot(FPulse*fs*2, fftshift(Exponential_Linear_Pulse_Spectrum),'--kp','LineWidth',1);
hold on
title('Frequency Response for ISI-Free Pulses with alpha='+string(alpha),'LineWidth',.1);
legend('RC', 'IPLCP', 'BTRC', 'ELP');
axis([-1.6 1.6 0 1])

grid on
xlabel('f/B')
ylabel('H(f)/T')
