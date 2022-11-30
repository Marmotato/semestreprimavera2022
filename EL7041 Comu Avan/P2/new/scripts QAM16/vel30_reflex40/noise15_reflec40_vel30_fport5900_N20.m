%Montecarlo 16QAM

%TO DO:

% % % % % % % % % % % % Initialization
clc;
clear all;

snrr=[-2:1:30];

noise = 15; %-5, 0 , 15, 30
reflec = 40; %5, 40
vel = 30; % 30, 80
fport = 5900; %700, 5900
N=20; %5, 10, 20

Bits = 10^5;
bps = 4;
%bits = randi([0 1],Bits,1);
%nb = Bits/2; %numero de símbolos, cada símbolo tamaño 2 
%nb=1593; %number of symbols symbols
nb = Bits/bps;

pilot=1+1i; % pilot symbol


%% 16QAM
% % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % %mapping 16QAM
inphase = [0.9 0.9 0.9 0.9 0.3 0.3 0.3 0.3];
quadr = [-0.9 -0.3 0.3 0.9 -0.9 -0.3 0.3 0.9];
inphase = [inphase;-inphase]; inphase = inphase(:);
quadr = [quadr;quadr]; quadr = quadr(:);
const = inphase + j*quadr;


M_16QAM = 16;
b_16QAM = randi([0 M_16QAM-1],nb,1); %generación random de símbolos
% % % % % % % % % % % % % % % % % % % %
tx=genqammod(b_16QAM,const); %general quadrature amplitude modulation

[Txp]=addpilot(tx,nb,pilot,N);%adding pilot and calculating no. of symbols per frame
len=length(Txp);
ct=rayleighfading(len, reflec, vel, fport); %genarating rayliegh fading channel

twn=Txp.*ct; %Multiplying Rayleigh channel coeeficients
% % % % % % % % % % % % % % noise generation
t=awgn(twn,noise,'measured','db' );%%%%% SNR
% % % % % % % % % % % % % % % % % % % % % % %
[rx rxp]=extractpilot(t,N,nb); %extracting pilot at receiver

%Channel Estimation CSI
csi=rxp/pilot; %calculating CSI of pilot
% % % % % % % % % % interpolation techniques
chnstinf=interpft(csi,nb);%%%%%%%% fft interpolation

t1 = 1:1:length(Txp);
t2 = 1:11:length(Txp); %Posición de las piloto
m = 1;
k = 0;
t3 = [];
for i = 1:11:length(Txp)
    t3 = [t3 t1(i+1:i+N)]; %Posición de los símbolos
    m = m + N;
    k = k + 1;
end
chnstinf2=interp1(t2,csi,t3,'spline');%%%%%%% spline cubic interpolation
chnstinf3=interp1(t2,csi,t3,'linear');%%%%%%% linear interpolation 
chnstinf4=interp1(t2,csi,t3,'pchip');%%%%%%%% cubic interpolation technique

% % % % % % % % % % %Calculating received signal --> ECUALIZACIÓN ZERO
% FORCING
for i=1:nb
RX(i)=rx(i)/chnstinf(i); %%%%calculating rxd symbols using fft interpolation technique
RX2(i)=rx(i)/chnstinf2(i);%%% using spline cubic interpolation
RX3(i)=rx(i)/chnstinf3(i);%%%%using linear interpolation
RX4(i)=rx(i)/chnstinf4(i);%%%% using cubic interpolation
end

[ch chp]=extpltcoef(ct,N,nb); %Extracting actual pilot coefficients and channel coefficients

%% PLOTS 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plots %%%%%%%%%%%%%%%%%%%%%%%

QAM_noise=awgn(Txp,noise,'measured','db' );%%%%% SNR
figure,plot(real(QAM_noise),imag(QAM_noise),'x');
%title('16QAM AFFECTED WITH NOISE (SNR='+string(noise)+')');
xlabel('REAL(DATA)');
ylabel('IMG(DATA)');

figure,plot(real(twn),imag(twn),'r.');
%title('16QAM AFFECTED WITH RAYLEIGH FADING');
xlabel('REAL(DATA)');
ylabel('IMG(DATA)');

figure,plot(real(rx),imag(rx),'r.');
%title('16QAM PLOT WITH NOISE AND RAYLEIGH FADING');
xlabel('REAL(DATA)');
ylabel('IMG(DATA)');

figure,
plot(real(RX),imag(RX),'r.');
hold on;
plot(real(tx),imag(tx),'b.');
grid on;
%title('16QAM PLOT');
xlabel('REAL(DATA)');
ylabel('IMG(DATA)');
legend('PLOT AT RX con Ecualización','PLOT AT TX' );
hold off

%% PLOTS 2
%%%%%%%%%%%%%plots 2%%%%%%%%%%%%%%%
%Función de transferencia del canal  y estimación del canal
cdb=10*log(abs(ch));
figure,plot(0:1:nb-1,cdb);
%title('Typical Rayleigh fading channel');
xlabel('t in samples');
ylabel('Power in dB');

figure,plot(0:1:nb-1,(abs(ch)),'r'); %Plotting actual value
hold on;
plot(0:1:nb-1,(abs(chnstinf)),'b'); %Plotting estimated value of fft
legend ('Actual value','Estimated value fft');
hold off;
%title('ACTUAL AND ESTIMATED CHANNEL COEFFICIENTS with FFT');
xlabel('Time in samples');
ylabel('Magnitude of coefficients');


figure,plot(0:1:nb-1,(abs(ch)),'r'); %Plotting actual value
hold on;
plot(0:1:nb-1,(abs(chnstinf2)),'b'); %Plotting estimated value by spline
legend ('Actual value','Estimated value spline');
hold off;
%title('ACTUAL AND ESTIMATED CHANNEL COEFFICIENTS with Spline');
xlabel('Time in samples');
ylabel('Magnitude of coefficients');

figure,plot(0:1:nb-1,(abs(ch)),'r'); %Plotting actual value
hold on;
plot(0:1:nb-1,(abs(chnstinf3)),'b'); %Plotting estimated value by linear
legend ('Actual value','Estimated value linear');
hold off;
%title('ACTUAL AND ESTIMATED CHANNEL COEFFICIENTS with Linear');
xlabel('Time in samples');
ylabel('Magnitude of coefficients');

figure,plot(0:1:nb-1,(abs(ch)),'r'); %Plotting actual value
hold on;
plot(0:1:nb-1,(abs(chnstinf4)),'b'); %Plotting estimated value by cubic
legend ('Actual value','Estimated value fft');
hold off;
%%title('ACTUAL AND ESTIMATED CHANNEL COEFFICIENTS with Cubic');
xlabel('Time in samples');
ylabel('Magnitude of coefficients');
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%% for BER curves%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% generating curves for snr ranging from -2 dB to 30dB

tot_rate1 = zeros(1,33);
tot_rate2 = zeros(1,33);
tot_rate3 = zeros(1,33);
tot_rate4 = zeros(1,33);

for asd=1:21
for l=1:length(snrr)
tt=awgn(twn,snrr(l),'measured','db' );
% % % % % % % % % % % % % % % % % % % % % % %
[rxt rxpt]=extractpilot(tt,N,nb); %extracting pilot at receiver
csit=rxpt/pilot; %calculating CSI of pilot
% % % % % % % % % % interpolation techniques

chnstinft=interpft(csit,nb);
t1t = 1:1:length(Txp);
t2t = 1:11:length(Txp); %Posición de las piloto
mt = 1;
kt = 0;
t3t = [];
for i = 1:11:length(Txp)
    t3t = [t3t t1(i+1:i+N)]; %Posición de los símbolos
    mt = mt + N;
    kt = kt + 1;
end
chnstinf2t=interp1(t2t,csit,t3t,'spline');
chnstinf3t=interp1(t2t,csit,t3t,'linear'); %DE AQUI SURGEN LOS NAN
chnstinf4t=interp1(t2t,csit,t3t,'pchip');
% % % % % % % % % % %Calculating received signal ECUALIZACIÓN ZF
for i=1:nb
RXt(i)=rxt(i)/chnstinft(i);
RX2t(i)=rxt(i)/chnstinf2t(i);
RX3t(i)=rxt(i)/chnstinf3t(i);
RX4t(i)=rxt(i)/chnstinf4t(i);
end


%Demodulación de la señal
rt=pskdemod(RXt,4);
r2t=pskdemod(RX2t,4);
r3t=pskdemod(RX3t,4);
r4t=pskdemod(RX4t,4);

% rt(isnan(rt))=0;
% r2t(isnan(r2t))=0;
 r3t(isnan(r3t))=0;
% r4t(isnan(r4t))=0;

[no_of_error1(l),rate1(l)]=biterr(b_16QAM.',rt) ; % error rate calculation for fft
[no_of_error2(l),rate2(l)]=biterr(b_16QAM.',r2t) ; % error rate calculation for spline
[no_of_error3(l),rate3(l)]=biterr(b_16QAM.',r3t) ; % error rate calculation for linear
[no_of_error4(l),rate4(l)]=biterr(b_16QAM.',r4t) ; % error rate calculation for cubic
end

tot_rate1 = tot_rate1 + rate1;
tot_rate2 = tot_rate2 + rate2;
tot_rate3 = tot_rate3 + rate3;
tot_rate4 = tot_rate4 + rate4;

end

tot_rate1 = tot_rate1/21;
tot_rate2 = tot_rate2/21;
tot_rate3 = tot_rate3/21;
tot_rate4 = tot_rate4/21;

% % % % % % % % % % % % % BER plot Monte Carlo % % % % % % % % % % %
figure
hold  on
semilogy(snrr, tot_rate1,'b-',snrr,tot_rate2,'r-',snrr,tot_rate3,'k-',snrr,tot_rate4,'g-');
legend('fft','cubic spline','linear','cubic');
%title('BER curves for different interpolation techniques');
xlabel('SNR in dB');
ylabel('BER');



% % % % % % % % % % % % % BER plot % % % % % % % % % % %
% figure
% hold  on
% semilogy(snrr, rate1,'b-',snrr,rate2,'r-',snrr,rate3,'k-',snrr,rate4,'g-');
% legend('fft','cubic spline','linear','cubic');
% %title('BER curves for different interpolation techniques');
% xlabel('SNR in dB');
% ylabel('BER');

% % % % % % % % % % % % % % % % % % % %
% % % % % % % 






%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%used functions
function [Txp]=addpilot(Tx,nsym,pilot,N) %Number of symbols per frame (Interpolation Factor)
m = 1;
k = 0;
for i = 1:N:nsym
    Txp(m) = pilot;
    Txp(m+1:m+N) = Tx(i:i+N-1);
    m = m + N + 1;
    k = k + 1;
end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove pilot coefficients -- We actually do not use this function.
function [ch chp]=extpltcoef(ct,N,nb)
m = 1;
ch = [];
for i=1:N+1:length(ct)
chp(m)=ct(i);
ch=[ch ct(i+1:i+N)];
m = m + 1;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Extraxting and removing pilot from transmitted signal
function [rx rxp]=extractpilot(t,N,nb)
m=1;
rx=[];
for i=1:N+1:length(t)
rxp(m)=t(i); %extracting only pilot symbols
rx=[rx t(i+1:i+N)]; %extracting data symbols
m=m+1;
end
end


%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Rayleigh fading
function [c]=rayleighfading(m, N, v, fport)
clc;
%N=5; % Number of reflections
c = 1.0793*10^9; %Velocidad de la luz en km/h
%v = 80; %Velocidad del móvil
f0 = fport*10^6; %Frecuencia portadora
fmax = v/c * f0;
A=1; %amplitude
t=0:1/f0:((m/10000)-(1/f0)); %sampling time
ct=zeros(1,m);
ph=2*pi* rand(1,N);
theta=2*pi*rand(1,N);
fd=fmax*cos(theta); %doppler shift
for n=1:m
for i=1:N
ct(n)=ct(n)+(A*exp(j*(2*pi*fd(i)*t(n)+ph(i))));

end
end
c=ct/sqrt(N); %channel coefficient
end