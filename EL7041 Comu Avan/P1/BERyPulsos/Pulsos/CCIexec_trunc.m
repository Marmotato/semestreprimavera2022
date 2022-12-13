%%execute all CCI

%% SIR = 10dB L=2
disp('alpha = 0.22, SIR = 10dB, L = 2, trunc = 5')

alpha = 0.22;
SIR = 10;
L = 2;
trunc = 5;

BER_CCI_i_truncated('RC', alpha, SIR, L, trunc)
BER_CCI_i_truncated('BTRC', alpha, SIR, L, trunc)
BER_CCI_i_truncated('IPLCP', alpha, SIR, L, trunc)
BER_CCI_i_truncated('ELP', alpha, SIR, L, trunc)

disp('alpha = 0.35, SIR = 10dB, L = 2')

alpha = 0.35;
SIR = 10;

BER_CCI_i_truncated('RC', alpha, SIR, L, trunc)
BER_CCI_i_truncated('BTRC', alpha, SIR, L, trunc)
BER_CCI_i_truncated('IPLCP', alpha, SIR, L, trunc)
BER_CCI_i_truncated('ELP', alpha, SIR, L, trunc)

alpha = 0.50;
SIR = 10;



disp('alpha = 0.5, SIR = 10dB, L = 2')

BER_CCI_i_truncated('RC', alpha, SIR, L, trunc)
BER_CCI_i_truncated('BTRC', alpha, SIR, L, trunc)
BER_CCI_i_truncated('IPLCP', alpha, SIR, L, trunc)
BER_CCI_i_truncated('ELP', alpha, SIR, L, trunc)

%% SIR = 10dB L=6

disp('alpha = 0.5, SIR = 10dB, L = 6')

alpha = 0.22;
SIR = 10;
L = 6;

BER_CCI_i_truncated('RC', alpha, SIR, L, trunc)
BER_CCI_i_truncated('BTRC', alpha, SIR, L, trunc)
BER_CCI_i_truncated('IPLCP', alpha, SIR, L, trunc)
BER_CCI_i_truncated('ELP', alpha, SIR, L, trunc)

disp('alpha = 0.35, SIR = 10dB, L = 6')

alpha = 0.35;
SIR = 10;

BER_CCI_i_truncated('RC', alpha, SIR, L, trunc)
BER_CCI_i_truncated('BTRC', alpha, SIR, L, trunc)
BER_CCI_i_truncated('IPLCP', alpha, SIR, L, trunc)
BER_CCI_i_truncated('ELP', alpha, SIR, L, trunc)

alpha = 0.50;
SIR = 10;



disp('alpha = 0.5, SIR = 10dB, L = 6')

BER_CCI_i_truncated('RC', alpha, SIR, L, trunc)
BER_CCI_i_truncated('BTRC', alpha, SIR, L, trunc)
BER_CCI_i_truncated('IPLCP', alpha, SIR, L, trunc)
BER_CCI_i_truncated('ELP', alpha, SIR, L, trunc)

%% SIR = 20dB L=2

disp('alpha = 0.22, SIR = 20dB, L = 2')

alpha = 0.22;
SIR = 20;

BER_CCI_i_truncated('RC', alpha, SIR, L, trunc)
BER_CCI_i_truncated('BTRC', alpha, SIR, L, trunc)
BER_CCI_i_truncated('IPLCP', alpha, SIR, L, trunc)
BER_CCI_i_truncated('ELP', alpha, SIR, L, trunc)


disp('alpha = 0.35, SIR = 20dB, L = 2')


alpha = 0.35;

BER_CCI_i_truncated('RC', alpha, SIR, L, trunc)
BER_CCI_i_truncated('BTRC', alpha, SIR, L, trunc)
BER_CCI_i_truncated('IPLCP', alpha, SIR, L, trunc)
BER_CCI_i_truncated('ELP', alpha, SIR, L, trunc)

alpha = 0.50;

disp('alpha = 0.5, SIR = 20dB, L = 2')

BER_CCI_i_truncated('RC', alpha, SIR, L, trunc)
BER_CCI_i_truncated('BTRC', alpha, SIR, L, trunc)
BER_CCI_i_truncated('IPLCP', alpha, SIR, L, trunc)
BER_CCI_i_truncated('ELP', alpha, SIR, L, trunc)

%% SIR = 20dB L=6

disp('alpha = 0.5, SIR = 20dB, L = 6')

SIR = 20;
alpha = 0.22;
L = 6;

BER_CCI_i_truncated('RC', alpha, SIR, L, trunc)
BER_CCI_i_truncated('BTRC', alpha, SIR, L, trunc)
BER_CCI_i_truncated('IPLCP', alpha, SIR, L, trunc)
BER_CCI_i_truncated('ELP', alpha, SIR, L, trunc)

disp('alpha = 0.35, SIR = 10dB, L = 6')

alpha = 0.35;

BER_CCI_i_truncated('RC', alpha, SIR, L, trunc)
BER_CCI_i_truncated('BTRC', alpha, SIR, L, trunc)
BER_CCI_i_truncated('IPLCP', alpha, SIR, L, trunc)
BER_CCI_i_truncated('ELP', alpha, SIR, L, trunc)

alpha = 0.50;

disp('alpha = 0.5, SIR = 20dB, L = 6')

BER_CCI_i_truncated('RC', alpha, SIR, L, trunc)
BER_CCI_i_truncated('BTRC', alpha, SIR, L, trunc)
BER_CCI_i_truncated('IPLCP', alpha, SIR, L, trunc)
BER_CCI_i_truncated('ELP', alpha, SIR, L, trunc)


%% SIR = 10dB L=2
disp('alpha = 0.22, SIR = 10dB, L = 2, trunc = 10')

alpha = 0.22;
SIR = 10;
L = 2;
trunc = 10;

BER_CCI_i_truncated('RC', alpha, SIR, L, trunc)
BER_CCI_i_truncated('BTRC', alpha, SIR, L, trunc)
BER_CCI_i_truncated('IPLCP', alpha, SIR, L, trunc)
BER_CCI_i_truncated('ELP', alpha, SIR, L, trunc)

disp('alpha = 0.35, SIR = 10dB, L = 2')

alpha = 0.35;
SIR = 10;

BER_CCI_i_truncated('RC', alpha, SIR, L, trunc)
BER_CCI_i_truncated('BTRC', alpha, SIR, L, trunc)
BER_CCI_i_truncated('IPLCP', alpha, SIR, L, trunc)
BER_CCI_i_truncated('ELP', alpha, SIR, L, trunc)

alpha = 0.50;
SIR = 10;



disp('alpha = 0.5, SIR = 10dB, L = 2')

BER_CCI_i_truncated('RC', alpha, SIR, L, trunc)
BER_CCI_i_truncated('BTRC', alpha, SIR, L, trunc)
BER_CCI_i_truncated('IPLCP', alpha, SIR, L, trunc)
BER_CCI_i_truncated('ELP', alpha, SIR, L, trunc)

%% SIR = 10dB L=6

disp('alpha = 0.5, SIR = 10dB, L = 6')

alpha = 0.22;
SIR = 10;
L = 6;

BER_CCI_i_truncated('RC', alpha, SIR, L, trunc)
BER_CCI_i_truncated('BTRC', alpha, SIR, L, trunc)
BER_CCI_i_truncated('IPLCP', alpha, SIR, L, trunc)
BER_CCI_i_truncated('ELP', alpha, SIR, L, trunc)

disp('alpha = 0.35, SIR = 10dB, L = 6')

alpha = 0.35;
SIR = 10;

BER_CCI_i_truncated('RC', alpha, SIR, L, trunc)
BER_CCI_i_truncated('BTRC', alpha, SIR, L, trunc)
BER_CCI_i_truncated('IPLCP', alpha, SIR, L, trunc)
BER_CCI_i_truncated('ELP', alpha, SIR, L, trunc)

alpha = 0.50;
SIR = 10;



disp('alpha = 0.5, SIR = 10dB, L = 6')

BER_CCI_i_truncated('RC', alpha, SIR, L, trunc)
BER_CCI_i_truncated('BTRC', alpha, SIR, L, trunc)
BER_CCI_i_truncated('IPLCP', alpha, SIR, L, trunc)
BER_CCI_i_truncated('ELP', alpha, SIR, L, trunc)

%% SIR = 20dB L=2

disp('alpha = 0.22, SIR = 20dB, L = 2')

alpha = 0.22;
SIR = 20;

BER_CCI_i_truncated('RC', alpha, SIR, L, trunc)
BER_CCI_i_truncated('BTRC', alpha, SIR, L, trunc)
BER_CCI_i_truncated('IPLCP', alpha, SIR, L, trunc)
BER_CCI_i_truncated('ELP', alpha, SIR, L, trunc)


disp('alpha = 0.35, SIR = 20dB, L = 2')


alpha = 0.35;

BER_CCI_i_truncated('RC', alpha, SIR, L, trunc)
BER_CCI_i_truncated('BTRC', alpha, SIR, L, trunc)
BER_CCI_i_truncated('IPLCP', alpha, SIR, L, trunc)
BER_CCI_i_truncated('ELP', alpha, SIR, L, trunc)

alpha = 0.50;

disp('alpha = 0.5, SIR = 20dB, L = 2')

BER_CCI_i_truncated('RC', alpha, SIR, L, trunc)
BER_CCI_i_truncated('BTRC', alpha, SIR, L, trunc)
BER_CCI_i_truncated('IPLCP', alpha, SIR, L, trunc)
BER_CCI_i_truncated('ELP', alpha, SIR, L, trunc)

%% SIR = 20dB L=6

disp('alpha = 0.5, SIR = 20dB, L = 6')

SIR = 20;
alpha = 0.22;
L = 6;

BER_CCI_i_truncated('RC', alpha, SIR, L, trunc)
BER_CCI_i_truncated('BTRC', alpha, SIR, L, trunc)
BER_CCI_i_truncated('IPLCP', alpha, SIR, L, trunc)
BER_CCI_i_truncated('ELP', alpha, SIR, L, trunc)

disp('alpha = 0.35, SIR = 10dB, L = 6')

alpha = 0.35;

BER_CCI_i_truncated('RC', alpha, SIR, L, trunc)
BER_CCI_i_truncated('BTRC', alpha, SIR, L, trunc)
BER_CCI_i_truncated('IPLCP', alpha, SIR, L, trunc)
BER_CCI_i_truncated('ELP', alpha, SIR, L, trunc)

alpha = 0.50;

disp('alpha = 0.5, SIR = 20dB, L = 6')

BER_CCI_i_truncated('RC', alpha, SIR, L, trunc)
BER_CCI_i_truncated('BTRC', alpha, SIR, L, trunc)
BER_CCI_i_truncated('IPLCP', alpha, SIR, L, trunc)
BER_CCI_i_truncated('ELP', alpha, SIR, L, trunc)
