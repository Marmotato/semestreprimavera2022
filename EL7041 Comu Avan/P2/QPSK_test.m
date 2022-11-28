N = 10^5; %número de bits
bits = randi([0 1],N,1);
SNRdB=[-2:1:30];

%% QPSK

M = 4; %quadrature

datos = randi([0 M-1], N, 1);%cantidad de datos modulados

tx = pskmod(datos, M, pi/M);
scatterplot(tx)
title('Constelación QPSK enviada');
xlabel('REAL(DATA)');
ylabel('IMG(DATA)');

for i=[-5, 0, 10, 30] %CONSTELACIONES CON EL EFECTO DE RUIDO
rx = awgn(tx, i); %Datos con ruido blanco agregado 
scatterplot(rx)
title('Constelación QPSK con ruido AWGN (SNR ='+string(i)+')');
xlabel('REAL(DATA)');
ylabel('IMG(DATA)');
end


