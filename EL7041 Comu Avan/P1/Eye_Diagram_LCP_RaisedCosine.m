
% Script for plotting the eye diagram of the system in which the transmitted signal is
% filtered by the RC and LCP.
% Added ELP and BTRC to filter.

clc
N  = 10^5; % number of symbols, solo se tomaran las primeras 1000000 muestras
am = 2*(rand(1,N)>0.5)-1 + (1i)*(2*(rand(1,N)>0.5)-1); %random binary sequence BPSK
fs = 10; % sampling frequency in Hz
alpha = 0.5; %roll-off factor
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%% RAISED COSINE FILTER

% defining the sinc filter
sincNum = sin(pi*[-fs:1/fs:fs]); % numerator of the sinc function
sincDen = (pi*[-fs:1/fs:fs]); % denominator of the sinc function
sincDenZero = find(abs(sincDen) < 10^-10);
sincOp = sincNum./sincDen;
sincOp(sincDenZero) = 1; % sin(pix/(pix) =1 for x =0

%Raised Cosine Filter
%alpha = 0.35;
cosNum = cos(alpha*pi*[-fs:1/fs:fs]);
cosDen = (1-(2*alpha*[-fs:1/fs:fs]).^2);
cosDenZero = find(abs(cosDen)<10^-10);
cosOp = cosNum./cosDen;
cosOp(cosDenZero) = pi/4;

gt_alpha35 = sincOp.*cosOp; %RC impulse responde
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

%% EXPONENTIAL LINEAR PULSE
constante=1; %alpha
beta = 0.1; %Beta


Exp2_2_Num=constante*sin(pi*[-fs:1/fs:fs]).*sin(pi*alpha*[-fs:1/fs:fs]); %PLP
Exp2_2_Den=pi^2*alpha*[-fs:1/fs:fs].^2;
Exp2_2DenZero=find(abs(Exp2_2_Den) < 10^-10);
Exp2_2_Total=Exp2_2_Num./Exp2_2_Den;
Exp2_Total(Exp2_2DenZero)=0.5;

Exp3 = exp(-pi*(beta/2)*[-fs:1/fs:fs].^2);

Exponential_Linear_Pulse=( Exp3.*Exp2_2_Total);
Exponential_Linear_Pulse(isnan(Exponential_Linear_Pulse))=1;

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


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%% upsampling the transmit sequence 
amUpSampled = [am;zeros(fs-1,length(am))];
amU = amUpSampled(:).';
%--------------------------------------------------------------------------
% filtered sequence using convolution fucntion
st_alpha35 = conv(amU,gt_alpha35);
IPLC=conv(amU,IPLC);
Exponential_Linear_Pulse = conv(amU,Exponential_Linear_Pulse);
BTRC=conv(amU,BTRC);
% %--------------------------------------------------------------------------
% %One only takes the first 1000000 samples for 
st_alpha35 = st_alpha35([1:1000000]);
IPLC=IPLC([1:1000000]);
Exponential_Linear_Pulse = Exponential_Linear_Pulse([1:1000000]);
BTRC = BTRC([1:1000000]);

st_alpha35_reshape = reshape(st_alpha35,fs*2,N*fs/20).';
IPLC_reshape=reshape(IPLC,fs*2,N*fs/20).';
Exponential_Linear_Pulse_reshape = reshape(Exponential_Linear_Pulse,fs*2,N*fs/20).';
BTRC_reshape=reshape(BTRC,fs*2,N*fs/20).';


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%% PLOT of Eye Diagrams
close all
figure(1);
plot([-0.99:1/fs:0.99],real(st_alpha35_reshape).','k');   
title('Eye Diagram of the RC pulse for an excess bandwidth alpha='+string(alpha),'LineWidth',.1);
xlabel('t/T')
ylabel('Amplitude') 
axis([-1 1 -2.5 2.5])
grid on
hold on


figure(2);
plot([-0.99:1/fs:0.99],real(IPLC_reshape).','k','LineWidth',.1);   
title('Eye Diagram of the IPLCP pulse for an excess bandwidth alpha='+string(alpha));
xlabel('t/T')
ylabel('Amplitude') 
axis([-1 1 -2.5 2.5])
grid on
hold on


figure(3);
plot([-0.99:1/fs:0.99],real(Exponential_Linear_Pulse_reshape).','k','LineWidth',.1);   
title('Eye Diagram of the ELP pulse for an excess bandwidth alpha='+string(alpha));
xlabel('t/T')
ylabel('Amplitude') 
axis([-1 1 -2.5 2.5])
grid on
hold on

figure(4);
plot([-0.99:1/fs:0.99],real(BTRC_reshape).','k','LineWidth',.1);   
title('Eye Diagram of the BTRC pulse for an excess bandwidth alpha='+string(alpha));
xlabel('t/T')
ylabel('Amplitude') 
axis([-1 1 -2.5 2.5])
grid on
hold on




