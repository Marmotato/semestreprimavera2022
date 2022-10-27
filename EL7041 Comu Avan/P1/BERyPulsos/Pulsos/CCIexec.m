%%execute all CCI

%% SIR = 10dB L=2
disp("alpha = 0.22, SIR = 10dB, L = 2")

alpha = 0.22;
SIR = 10;
L = 2;

BER_CCI_i('RC', alpha, SIR, L)
BER_CCI_i('BTRC', alpha, SIR, L)
BER_CCI_i('IPLCP', alpha, SIR, L)
BER_CCI_i('ELP', alpha, SIR, L)

disp("alpha = 0.35, SIR = 10dB, L = 2")

alpha = 0.35;
SIR = 10;

BER_CCI_i('RC', alpha, SIR, L)
BER_CCI_i('BTRC', alpha, SIR, L)
BER_CCI_i('IPLCP', alpha, SIR, L)
BER_CCI_i('ELP', alpha, SIR, L)

alpha = 0.50;
SIR = 10;



disp("alpha = 0.5, SIR = 10dB, L = 2")

BER_CCI_i('RC', alpha, SIR, L)
BER_CCI_i('BTRC', alpha, SIR, L)
BER_CCI_i('IPLCP', alpha, SIR, L)
BER_CCI_i('ELP', alpha, SIR, L)

%% SIR = 10dB L=6

disp("alpha = 0.5, SIR = 10dB, L = 6")

alpha = 0.22;
SIR = 10;
L = 6;

BER_CCI_i('RC', alpha, SIR, L)
BER_CCI_i('BTRC', alpha, SIR, L)
BER_CCI_i('IPLCP', alpha, SIR, L)
BER_CCI_i('ELP', alpha, SIR, L)

disp("alpha = 0.35, SIR = 10dB, L = 6")

alpha = 0.35;
SIR = 10;

BER_CCI_i('RC', alpha, SIR, L)
BER_CCI_i('BTRC', alpha, SIR, L)
BER_CCI_i('IPLCP', alpha, SIR, L)
BER_CCI_i('ELP', alpha, SIR, L)

alpha = 0.50;
SIR = 10;



disp("alpha = 0.5, SIR = 10dB, L = 6")

BER_CCI_i('RC', alpha, SIR, L)
BER_CCI_i('BTRC', alpha, SIR, L)
BER_CCI_i('IPLCP', alpha, SIR, L)
BER_CCI_i('ELP', alpha, SIR, L)

%% SIR = 20dB L=2

disp("alpha = 0.22, SIR = 20dB, L = 2")

alpha = 0.22;
SIR = 20;

BER_CCI_i('RC', alpha, SIR, L)
BER_CCI_i('BTRC', alpha, SIR, L)
BER_CCI_i('IPLCP', alpha, SIR, L)
BER_CCI_i('ELP', alpha, SIR, L)


disp("alpha = 0.35, SIR = 20dB, L = 2")


alpha = 0.35;

BER_CCI_i('RC', alpha, SIR, L)
BER_CCI_i('BTRC', alpha, SIR, L)
BER_CCI_i('IPLCP', alpha, SIR, L)
BER_CCI_i('ELP', alpha, SIR, L)

alpha = 0.50;

disp("alpha = 0.5, SIR = 20dB, L = 2")

BER_CCI_i('RC', alpha, SIR, L)
BER_CCI_i('BTRC', alpha, SIR, L)
BER_CCI_i('IPLCP', alpha, SIR, L)
BER_CCI_i('ELP', alpha, SIR, L)

%% SIR = 20dB L=6

disp("alpha = 0.5, SIR = 20dB, L = 6")

SIR = 20
alpha = 0.22;
L = 6;

BER_CCI_i('RC', alpha, SIR, L)
BER_CCI_i('BTRC', alpha, SIR, L)
BER_CCI_i('IPLCP', alpha, SIR, L)
BER_CCI_i('ELP', alpha, SIR, L)

disp("alpha = 0.35, SIR = 10dB, L = 6")

alpha = 0.35;

BER_CCI_i('RC', alpha, SIR, L)
BER_CCI_i('BTRC', alpha, SIR, L)
BER_CCI_i('IPLCP', alpha, SIR, L)
BER_CCI_i('ELP', alpha, SIR, L)

alpha = 0.50;

disp("alpha = 0.5, SIR = 20dB, L = 6")

BER_CCI_i('RC', alpha, SIR, L)
BER_CCI_i('BTRC', alpha, SIR, L)
BER_CCI_i('IPLCP', alpha, SIR, L)
BER_CCI_i('ELP', alpha, SIR, L)

