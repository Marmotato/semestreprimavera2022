
% Script for plotting the eye diagram of the system in which the transmitted signal is
% filtered by the RC and LCP.
% Added ELP and BTRC to filter.

clc
N  = 10^5; % number of symbols, solo se tomaran las primeras 1000000 muestras
am = 2*(rand(1,N)>0.5)-1 + (1i)*(2*(rand(1,N)>0.5)-1); %random binary sequence BPSK
fs = 10; % sampling frequency in Hz
alpha = 0.22; %roll-off factor
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%RAISED COSINE FILTER

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
%Linear Combination Pulse
% constante*PLP(n=1) +(1-constante)RC
constante=1.7;  %valor de la constante Beta --> referencia 

Exp1_Num=((1-constante)*sin(pi*[-fs:1/fs:fs]).*cos(pi*alpha*[-fs:1/fs:fs])); %RC
Exp1_Den=pi*[-fs:1/fs:fs].*(1-(2*alpha*[-fs:1/fs:fs]).^2);
Exp1DenZero=find(abs(Exp1_Den) < 10^-10);
Exp1_Total=Exp1_Num./Exp1_Den;
Exp1_Total(Exp1DenZero)=0;

Exp2_Num=constante*sin(pi*[-fs:1/fs:fs]).*sin(pi*alpha*[-fs:1/fs:fs]); %PLP
Exp2_Den=pi^2*alpha*[-fs:1/fs:fs].^2;
Exp2DenZero=find(abs(Exp2_Den) < 10^-10);
Exp2_Total=Exp2_Num./Exp2_Den;
Exp2_Total(Exp2DenZero)=1;

Linear_Combination_Pulse_PLP=(Exp1_Total+Exp2_Total);

%EXPONENTIAL LINEAR PULSE
constante=1; %alpha
beta = 0.1; %Beta


Exp2_2_Num=constante*sin(pi*[-fs:1/fs:fs]).*sin(pi*alpha*[-fs:1/fs:fs]); %PLP
Exp2_2_Den=pi^2*alpha*[-fs:1/fs:fs].^2;
Exp2_2DenZero=find(abs(Exp2_Den) < 10^-10);
Exp2_2_Total=Exp2_2_Num./Exp2_Den;
Exp2_Total(Exp2_2DenZero)=0.5;

Exp3 = exp(-pi*(beta/2)*[-fs:1/fs:fs].^2);

Exponential_Linear_Pulse=( Exp3.*Exp2_2_Total);
Exponential_Linear_Pulse(isnan(Exponential_Linear_Pulse))=1;

%BETTER THAN RAISED COSINE PULSE

beta2 = (pi * alpha) / log(2);

% Calculate the sinc expression
sinc_BTRC_num = sin([-fs:1/fs:fs]);
sinc_BTRC_den = [-fs:1/fs:fs];
sinc_BTRC_tot = sinc_BTRC_num./sinc_BTRC_den;
sinc_BTRC_tot_zero=find(abs(sinc_BTRC_tot) < 10^-10);

% Calculate the second expression for BTRC
numerador_second_BTRC = (2 * beta2 * [-fs:1/fs:fs]) .* sin(pi * alpha * [-fs:1/fs:fs]);
numerador_second_BTRC = numerador_second_BTRC + (2 .* cos(pi * alpha * [-fs:1/fs:fs]) - 1);

% Never 0
denominador_second_BTRC = (beta2 * [-fs:1/fs:fs]).^2 + 1;

% Get the BTRC Pulse
BTRC = ((sinc_BTRC_tot .* numerador_second_BTRC) ./ denominador_second_BTRC);
%BTRC = arrayfun(@(x) BTRC(x,alpha),[-fs:1/fs:fs]);
BTRC(sinc_BTRC_tot_zero) = (numerador_second_BTRC(sinc_BTRC_tot_zero) / denominador_second_BTRC(sinc_BTRC_tot_zero));

BTRC(isnan(BTRC))=1;


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% upsampling the transmit sequence 
amUpSampled = [am;zeros(fs-1,length(am))];
amU = amUpSampled(:).';
%--------------------------------------------------------------------------
% filtered sequence using convolution fucntion
st_alpha35 = conv(amU,gt_alpha35);
Linear_Combination_Pulse_PLP=conv(amU,Linear_Combination_Pulse_PLP);
Exponential_Linear_Pulse = conv(amU,Exponential_Linear_Pulse);
BTRC=conv(amU,BTRC);
% %--------------------------------------------------------------------------
% %One only takes the first 1000000 samples for 
st_alpha35 = st_alpha35([1:1000000]);
Linear_Combination_Pulse_PLP=Linear_Combination_Pulse_PLP([1:1000000]);
Exponential_Linear_Pulse = Exponential_Linear_Pulse([1:1000000]);
BTRC = BTRC([1:1000000]);

st_alpha35_reshape = reshape(st_alpha35,fs*2,N*fs/20).';
Linear_Combination_Pulse_PLP_reshape=reshape(Linear_Combination_Pulse_PLP,fs*2,N*fs/20).';
Exponential_Linear_Pulse_reshape = reshape(Exponential_Linear_Pulse,fs*2,N*fs/20).';
BTRC_reshape=reshape(BTRC,fs*2,N*fs/20).';


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%PLOT of Eye Diagrams
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
plot([-0.99:1/fs:0.99],real(Linear_Combination_Pulse_PLP_reshape).','k','LineWidth',.1);   
title('Eye Diagram of the LCP pulse for an excess bandwidth alpha='+string(alpha));
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




