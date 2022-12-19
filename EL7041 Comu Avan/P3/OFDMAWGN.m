%% Se inicializa OFDM

%%Parámetros de simulación para OFDM BPSK, enviando 10^5 símbolos OFDM

reflex_l = [40];
fcs_l=[700, 5900];
Ns_l = [5,10,20];

for i1 = 1:length(reflex_l)
for j1 = 1:length(fcs_l)
for p1 = 1:3

nFFT = 64; %Tamaño de la FFT
nDSC = 52; %Número de portadoras con información
nBitPerSym = 52; %Número de símbolos por portadora (BPSK 1 bit por símbolo)
nSym = 10^5; %Número de símbolos transmitidos en OFDM
SNR = -5; %Valor de SNR en dB
N = Ns_l(p1); %Número de símbolos por piloto
pilot = 1+1i; %Señal piloto
reflx = reflex_l(i1); %Número de reflexiones
vel = 80; %Velocidad en km/h
fc = fcs_l(j1); %Frecuencia de la portadora en MHz


%Modulación

ipBit = randi([0 1], 1, nBitPerSym*nSym); %Ceros y unos aleatorios
ipMod = 2*ipBit - 1; %Se realiza la modulación BPSK
ipMod = reshape(ipMod, nBitPerSym, nSym); %Eje vertical: frecuencias; eje 
%horizontal: símbolos

%Se asignan los ceros a las bandas de guarda y la DC subcarrier
xF = [zeros(6, nSym); ipMod(1:nBitPerSym/2, :); zeros(1, nSym);...
    ipMod(nBitPerSym/2 + 1:nBitPerSym, :); zeros(5, nSym)];

%Se toma la transformada inversa de Fourier para llevar los símbolos a
%banda base en un DeltaF (exp(1i*2*pi*DeltaF*t)), normalizada
xt = (nFFT/sqrt(nDSC))*ifft(xF, [], 1);

%Se añade el prefijo cíclico
xt = [xt(49:64, :); xt];

%Se añade la señal piloto

TxPilot = zeros(80, round((N + 1)/N * nSym), 1); %Reserva espacio en memoria

m = 1;
for j = 1:N:nSym
    TxPilot(:, m) = pilot;
    TxPilot(:, m + 1:m + N) = xt(:, j:j + N - 1);
    m = m + N + 1;
end

%Se convierte xt en una sola trama de símbolos concatenándolos
TxPilot = reshape(TxPilot, 1, []);

%fig = figure('Position', get(0, 'ScreenSize'));
%plot(real(ipMod), imag(ipMod), 'x')
%xlabel('Componente en fase')
%ylabel('Componente en cuadratura')
%grid on
%grid minor
%fontsize(fig, 25, 'pixels')
%set(fig, 'color', 'w');
%saveas(fig, 'Constelación transmitida.png', 'png')

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

chRayleighTotal = sum(chRayleighTotal)/sqrt(reflx); %Respuesta del canal
TxRayleigh = TxPilot.*chRayleighTotal; %Señal transmitida por el canal

RxRayleighTotal = awgn(TxRayleigh, SNR, 'measured', 'dB'); %Se añade AWGN
RxPilotRayleigh = zeros(1, numel(xt)/N); %Reserva espacio en memoria
RxRayleigh = zeros(1, numel(xt)); %Reserva espacio en memoria
chResponseRayleighPilot = zeros(1, numel(xt)/N); %Reserva espacio en memoria
chResponseRayleigh = zeros(1, numel(xt)); %Reserva espacio en memoria

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
t3 = zeros(numel(xt), 1); %Reserva espacio en memoria
m = 1;
for i = 1:N + 1:length(RxRayleighTotal)
    t3(m:m + N - 1) = t1(i + 1:i + N);
    m = m + N;
end

chinf1 = interpft(CSI, numel(xt)); %Interpolación FFT
chinf2 = interp1(t2, CSI, t3, 'spline').'; %Interpolación spline
chinf3 = interp1(t2, CSI, t3, 'linear').'; %Interpolación lineal
chinf4 = interp1(t2, CSI, t3, 'cubic').'; %Interpolación cúbica

%Se hace una ecualización Zero-Forcing

RX1 = RxRayleigh./chinf1;
RX2 = RxRayleigh./chinf2;
RX3 = RxRayleigh./chinf3;
RX4 = RxRayleigh./chinf4;


%% Receptor

%Se devuelve la señal a símbolos
RX1 = reshape(RX1, 80, []);
RX2 = reshape(RX2, 80, []);
RX3 = reshape(RX3, 80, []);
RX4 = reshape(RX4, 80, []);

%Se le quitan los prefijos cíclicos
RX1 = RX1(17:end, :);
RX2 = RX2(17:end, :);
RX3 = RX3(17:end, :);
RX4 = RX4(17:end, :);

%Se convierte la señal recibida al dominio de la frecuencia nuevamente
RF1 = (sqrt(nDSC)/nFFT)*fftshift(fft(RX1, [], 1));
RF2 = (sqrt(nDSC)/nFFT)*fftshift(fft(RX2, [], 1));
RF3 = (sqrt(nDSC)/nFFT)*fftshift(fft(RX3, [], 1));
RF4 = (sqrt(nDSC)/nFFT)*fftshift(fft(RX4, [], 1));

%Se quitan las bandas de guarda
RMod1 = RF1([7:6+nBitPerSym/2 7+nBitPerSym/2+1:7+nBitPerSym], :);
RMod2 = RF2([7:6+nBitPerSym/2 7+nBitPerSym/2+1:7+nBitPerSym], :);
RMod3 = RF3([7:6+nBitPerSym/2 7+nBitPerSym/2+1:7+nBitPerSym], :);
RMod4 = RF4([7:6+nBitPerSym/2 7+nBitPerSym/2+1:7+nBitPerSym], :);

%Demodulación BPSK
ipModHat1 = 2*floor(real(RMod1/2)) + 1;
ipModHat1(find(ipModHat1 > 1)) = +1;
ipModHat1(find(ipModHat1 < -1)) = -1;

ipModHat2 = 2*floor(real(RMod2/2)) + 1;
ipModHat2(find(ipModHat2 > 1)) = +1;
ipModHat2(find(ipModHat2 < -1)) = -1;

ipModHat3 = 2*floor(real(RMod3/2)) + 1;
ipModHat3(find(ipModHat3 > 1)) = +1;
ipModHat3(find(ipModHat3 < -1)) = -1;

ipModHat4 = 2*floor(real(RMod4/2)) + 1;
ipModHat4(find(ipModHat4 > 1)) = +1;
ipModHat4(find(ipModHat4 < -1)) = -1;

%Se convierten los valores modulados en analógico a bits
ipBitHat1 = (ipModHat1 + 1)/2;
ipBitHat1 = reshape(ipBitHat1, 1, []);
ipBitHat1(isnan(ipBitHat1)) = 0;

ipBitHat2 = (ipModHat2 + 1)/2;
ipBitHat2 = reshape(ipBitHat2, 1, []);
ipBitHat2(isnan(ipBitHat2)) = 0;

ipBitHat3 = (ipModHat3 + 1)/2;
ipBitHat3 = reshape(ipBitHat3, 1, []);
ipBitHat3(isnan(ipBitHat3)) = 0;

ipBitHat4 = (ipModHat4 + 1)/2;
ipBitHat4 = reshape(ipBitHat4, 1, []);
ipBitHat4(isnan(ipBitHat4)) = 0;

%Se cuenta el número de errores
nErr1 = biterr(ipBitHat1, ipBit);
nErr2 = biterr(ipBitHat2, ipBit);
nErr3 = biterr(ipBitHat3, ipBit);
nErr4 = biterr(ipBitHat4, ipBit);

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

for k = 1:1
    for l = 1:length(SNR_dB)
        RxRayleighTotal = awgn(TxRayleigh, SNR_dB(l), 'measured', 'dB'); %Se añade AWGN
        RxPilotRayleigh = zeros(1, numel(xt)/N); %Reserva espacio en memoria
        RxRayleigh = zeros(1, numel(xt)); %Reserva espacio en memoria
        chResponseRayleighPilot = zeros(1, numel(xt)/N); %Reserva espacio en memoria
        chResponseRayleigh = zeros(1, numel(xt)); %Reserva espacio en memoria

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
        t3 = zeros(numel(xt), 1); %Reserva espacio en memoria
        m = 1;
        for i = 1:N + 1:length(RxRayleighTotal)
            t3(m:m + N - 1) = t1(i + 1:i + N);
            m = m + N;
        end

        chinf1 = interpft(CSI, numel(xt)); %Interpolación FFT
        chinf2 = interp1(t2, CSI, t3, 'spline').'; %Interpolación spline
        chinf3 = interp1(t2, CSI, t3, 'linear').'; %Interpolación lineal
        chinf4 = interp1(t2, CSI, t3, 'cubic').'; %Interpolación cúbica

        %Se hace una ecualización Zero-Forcing

        RX1 = RxRayleigh./chinf1;
        RX2 = RxRayleigh./chinf2;
        RX3 = RxRayleigh./chinf3;
        RX4 = RxRayleigh./chinf4;

        %Receptor

        %Se devuelve la señal a símbolos
        RX1 = reshape(RX1, 80, []);
        RX2 = reshape(RX2, 80, []);
        RX3 = reshape(RX3, 80, []);
        RX4 = reshape(RX4, 80, []);

        %Se le quitan los prefijos cíclicos
        RX1 = RX1(17:end, :);
        RX2 = RX2(17:end, :);
        RX3 = RX3(17:end, :);
        RX4 = RX4(17:end, :);

        %Se convierte la señal recibida al dominio de la frecuencia nuevamente
        RF1 = (sqrt(nDSC)/nFFT)*fftshift(fft(RX1, [], 1));
        RF2 = (sqrt(nDSC)/nFFT)*fftshift(fft(RX2, [], 1));
        RF3 = (sqrt(nDSC)/nFFT)*fftshift(fft(RX3, [], 1));
        RF4 = (sqrt(nDSC)/nFFT)*fftshift(fft(RX4, [], 1));

        %Se quitan las bandas de guarda
        RMod1 = RF1([7:6+nBitPerSym/2 7+nBitPerSym/2+1:7+nBitPerSym], :);
        RMod2 = RF2([7:6+nBitPerSym/2 7+nBitPerSym/2+1:7+nBitPerSym], :);
        RMod3 = RF3([7:6+nBitPerSym/2 7+nBitPerSym/2+1:7+nBitPerSym], :);
        RMod4 = RF4([7:6+nBitPerSym/2 7+nBitPerSym/2+1:7+nBitPerSym], :);

        %Demodulación BPSK
        ipModHat1 = 2*floor(real(RMod1/2)) + 1;
        ipModHat1(ipModHat1 > 1) = +1;
        ipModHat1(ipModHat1 < -1) = -1;

        ipModHat2 = 2*floor(real(RMod2/2)) + 1;
        ipModHat2(ipModHat2 > 1) = +1;
        ipModHat2(ipModHat2 < -1) = -1;

        ipModHat3 = 2*floor(real(RMod3/2)) + 1;
        ipModHat3(ipModHat3 > 1) = +1;
        ipModHat3(ipModHat3 < -1) = -1;

        ipModHat4 = 2*floor(real(RMod4/2)) + 1;
        ipModHat4(ipModHat4 > 1) = +1;
        ipModHat4(ipModHat4 < -1) = -1;

        %Se convierten los valores modulados en analógico a bits
        ipBitHat1 = (ipModHat1 + 1)/2;
        ipBitHat1 = reshape(ipBitHat1, 1, []);
        ipBitHat1(isnan(ipBitHat1)) = 0;

        ipBitHat2 = (ipModHat2 + 1)/2;
        ipBitHat2 = reshape(ipBitHat2, 1, []);
        ipBitHat2(isnan(ipBitHat2)) = 0;

        ipBitHat3 = (ipModHat3 + 1)/2;
        ipBitHat3 = reshape(ipBitHat3, 1, []);
        ipBitHat3(isnan(ipBitHat3)) = 0;

        ipBitHat4 = (ipModHat4 + 1)/2;
        ipBitHat4 = reshape(ipBitHat4, 1, []);
        ipBitHat4(isnan(ipBitHat4)) = 0;
        
        %Se calculan los BER
    
        [nerror1(l), prob1(l)] = biterr(ipBit, ipBitHat1); %BER para FFT
        [nerror2(l), prob2(l)] = biterr(ipBit, ipBitHat2); %BER para Spline
        [nerror3(l), prob3(l)] = biterr(ipBit, ipBitHat3); %BER para Lineal
        [nerror4(l), prob4(l)] = biterr(ipBit, ipBitHat4); %BER para Cúbica
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

fig = figure('Position', get(0, 'ScreenSize'));
semilogy(SNR_dB, BERFFT,SNR_dB, BERSpline, SNR_dB, BERLineal, SNR_dB, BERCubic);
legend('FFT','Spline','Lineal','Cúbica');
xlabel('SNR (dB)')
ylabel('BER')
grid on
grid minor
set(fig, 'color', 'w');
saveas(fig, 'BER N = ' + string(N) + ', ' +...
    string(reflx) + ' reflexiones, ' + string(vel) + ' kmh, ' + string(fc)...
    + ' MHz.png', 'png')

end
end
end