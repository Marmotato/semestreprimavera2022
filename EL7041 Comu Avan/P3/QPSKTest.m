%Se quiere calcular el BER de QPSK, enviando 2^20 bits

nbits = 100000; %Bits a enviar
bps = 2; %Bits por símbolo
nsymb = round(nbits/bps); %Número total de símbolos a enviar
pilot = 1 + 1i; %Señal piloto

%Parámetros de simulación

SNR = 30; %SNR en dB
N = 5; %Cantidad de símbolos por frame
reflx = 40; %Número de reflexiones
vel = 80; %Velocidad relativa en km/h
fc = 5900*10^6; %Frecuencia de la portadora central en Hz

%Se genera la señal a transmitir

M = 2^bps; %Grado de modulación
bits = randi([0 M-1], nsymb, 1); %Bits aleatorios a enviar

qpskmod = comm.QPSKModulator; %Objeto modulador
qpskdemod = comm.QPSKDemodulator; %Objero demodulador
Tx = qpskmod(bits); %Señal a transmitir

%Se añaden las señales piloto

TxPilot = zeros(round((N + 1)/N * nsymb), 1); %Reserva espacio en memoria

m = 1;
for i = 1:N:nsymb
    TxPilot(m) = pilot;
    TxPilot(m + 1:m + N) = Tx(i:i + N - 1);
    m = m + N + 1;
end

%% Canal Rayleigh

A = 1; %Ganancia de cada path
L = length(TxPilot); %Tamaño de la muestra
fs = 100000; %Frecuencia de muestreo
t = linspace(0, L/fs, L); %Tiempo de respuesta
chRayleighTotal = zeros(N, L); %Reserva espacio en memoria
phase = 2*pi*rand(1, reflx); %Fase de las reflexiones
theta = 2*pi*rand(1, reflx); %Ángulo del Doppler
c = 1.0793*10^9; %Velocidad de la luz en km/h
fd = vel/c * fc * cos(theta); %Desviación por doppler
for i = 1:L
    for j = 1:reflx
        chRayleighTotal(j, i) = A*exp(1i*(2*pi*fd(j)*t(i) + phase(j)));
    end
end

chRayleighTotal = sum(chRayleighTotal)'/sqrt(reflx); %Respuesta del canal
TxRayleigh = TxPilot.*chRayleighTotal; %Señal transmitida por el canal

RxRayleighTotal = awgn(TxRayleigh, SNR, 'measured', 'dB'); %Se añade AWGN
RxPilotRayleigh = zeros(round(nsymb/N), 1); %Reserva espacio en memoria
RxRayleigh = zeros(round(nsymb), 1); %Reserva espacio en memoria
chResponseRayleighPilot = zeros(round(nsymb/N), 1); %Reserva espacio en memoria
chResponseRayleigh = zeros(nsymb, 1); %Reserva espacio en memoria

m = 1;
n = 1;
for i = 1:N + 1:length(RxRayleighTotal)
    RxPilotRayleigh(m) = RxRayleighTotal(i);
    RxRayleigh(n:n + N - 1) = RxRayleighTotal(i + 1:i + N);
    chResponseRayleighPilot(m) = chRayleighTotal(i);
    chResponseRayleigh(n:n + N - 1) = chRayleighTotal(i + 1:i + N);
    m = m + 1;
    n = n + N;
end

%Estimación perfecta de canal Rayleigh

CSI = RxPilotRayleigh/pilot;

t1 = 1:length(RxRayleighTotal); %Largo total de las señales
t2 = 1:N + 1:length(RxRayleighTotal); %Posición de las señales piloto
t3 = zeros(nsymb, 1); %Reserva espacio en memoria
m = 1;
for i = 1:N + 1:length(RxRayleighTotal)
    t3(m:m + N - 1) = t1(i + 1:i + N);
    m = m + N;
end

chinf1 = interpft(CSI, nsymb); %Interpolación FFT
chinf2 = interp1(t2, CSI, t3, 'spline'); %Interpolación spline
chinf3 = interp1(t2, CSI, t3, 'linear'); %Interpolación lineal
chinf4 = interp1(t2, CSI, t3, 'cubic'); %Interpolación cúbica

%Se hace una ecualización Zero-Forcing

RX1 = RxRayleigh./chinf1;
RX2 = RxRayleigh./chinf2;
RX3 = RxRayleigh./chinf3;
RX4 = RxRayleigh./chinf4;

%% Figuras

%Señal recibida en canal AWGN
RxAWGN = awgn(TxPilot, SNR, 'measured', 'dB');
figure()
plot(RxAWGN, 'x')

%Señal recibida sólo canal Rayleigh
figure()
plot(TxRayleigh, 'x');

%Señal recibida en canal rayleigh con AWGN
figure()
plot(RxRayleigh, 'x')

%Señal recibida con ecualización ZF
figure()
plot(RX1, 'r.')
hold on
plot(Tx, 'o')

%Función de transferencia del canal
chRayleighdB = 20*log10(abs(chResponseRayleigh));
figure()
plot(chRayleighdB)

%Estimación FFT
figure()
plot(abs(chResponseRayleigh), 'r')
hold on
plot(abs(chinf1), 'b')

%Estimación Spline
figure()
plot(abs(chResponseRayleigh), 'r')
hold on
plot(abs(chinf2), 'b')

%Estimación Lineal
figure()
plot(abs(chResponseRayleigh), 'r')
hold on
plot(abs(chinf3), 'b')

%Estimación Cúbica
figure()
plot(abs(chResponseRayleigh), 'r')
hold on
plot(abs(chinf4), 'b')

%% Curvas de BER en AWGN

SNR_dB = -2:30; %SNR en dB

nerror1 = zeros(1, length(SNR_dB));
nerror2 = zeros(1, length(SNR_dB));
nerror3 = zeros(1, length(SNR_dB));
nerror4 = zeros(1, length(SNR_dB));
prob1 = zeros(1, length(SNR_dB));
prob2 = zeros(1, length(SNR_dB));
prob3 = zeros(1, length(SNR_dB));
prob4 = zeros(1, length(SNR_dB));

%Montecarlo

BERFFT = zeros(21, 33);
BERSpline = zeros(21, 33);
BERLineal = zeros(21, 33);
BERCubic = zeros(21, 33);

for k = 1:21
    for l = 1:length(SNR_dB)
        RxRayleighTotal = awgn(TxRayleigh, SNR_dB(l), 'measured', 'dB'); %Se añade AWGN
        RxPilotRayleigh = zeros(round(nsymb/N), 1); %Reserva espacio en memoria
        RxRayleigh = zeros(round(nsymb), 1); %Reserva espacio en memoria
        chResponseRayleighPilot = zeros(round(nsymb/N), 1); %Reserva espacio en memoria
        chResponseRayleigh = zeros(nsymb, 1); %Reserva espacio en memoria

        m = 1;
        n = 1;
        for i = 1:N + 1:length(RxRayleighTotal)
            RxPilotRayleigh(m) = RxRayleighTotal(i);
            RxRayleigh(n:n + N - 1) = RxRayleighTotal(i + 1:i + N);
            chResponseRayleighPilot(m) = chRayleighTotal(i);
            chResponseRayleigh(n:n + N - 1) = chRayleighTotal(i + 1:i + N);
            m = m + 1;
            n = n + N;
        end

        %Estimación perfecta de canal Rayleigh

        CSI = RxPilotRayleigh/pilot;

        t1 = 1:length(RxRayleighTotal); %Largo total de las señales
        t2 = 1:N + 1:length(RxRayleighTotal); %Posición de las señales piloto
        t3 = zeros(nsymb, 1); %Reserva espacio en memoria
        m = 1;
        for i = 1:N + 1:length(RxRayleighTotal)
            t3(m:m + N - 1) = t1(i + 1:i + N);
            m = m + N;
        end

        chinf1 = interpft(CSI, nsymb); %Interpolación FFT
        chinf2 = interp1(t2, CSI, t3, 'spline'); %Interpolación spline
        chinf3 = interp1(t2, CSI, t3, 'linear'); %Interpolación lineal
        chinf4 = interp1(t2, CSI, t3, 'cubic'); %Interpolación cúbica

        %Se hace una ecualización Zero-Forcing

        RX1 = RxRayleigh./chinf1;
        RX2 = RxRayleigh./chinf2;
        RX3 = RxRayleigh./chinf3;
        RX4 = RxRayleigh./chinf4;

        %Se cambian los nan por 0

        RX1(isnan(RX1)) = 0;
        RX2(isnan(RX2)) = 0;
        RX3(isnan(RX3)) = 0;
        RX4(isnan(RX4)) = 0;

        %Se demodula la señal

        r1 = qpskdemod(RX1);
        r2 = qpskdemod(RX2);
        r3 = qpskdemod(RX3);
        r4 = qpskdemod(RX4);

        %Se cambian los nan por 0

        r1(isnan(r1)) = 0;
        r2(isnan(r2)) = 0;
        r3(isnan(r3)) = 0;
        r4(isnan(r4)) = 0;

        %Se calculan los BER
    
        [nerror1(l), prob1(l)] = biterr(bits, r1); %BER para FFT
        [nerror2(l), prob2(l)] = biterr(bits, r2); %BER para Spline
        [nerror3(l), prob3(l)] = biterr(bits, r3); %BER para Lineal
        [nerror4(l), prob4(l)] = biterr(bits, r4); %BER para Cúbica
    end
    BERFFT(k, :) = prob1;
    BERSpline(k, :) = prob2;
    BERLineal(k, :) = prob3;
    BERCubic(k, :) = prob4;
end

BERFFT = mean(BERFFT);
BERSpline = mean(BERSpline);
BERLineal = mean(BERLineal);
BERCubic = mean(BERCubic);

figure()
semilogy(SNR_dB, BERFFT,SNR_dB, BERSpline, SNR_dB, BERLineal, SNR_dB, BERCubic);
legend('FFT','Spline','Lineal','Cúbica');
