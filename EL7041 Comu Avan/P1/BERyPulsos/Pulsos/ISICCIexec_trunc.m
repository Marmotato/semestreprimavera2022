%%execute all CCI and ISI

%% SIR = 15dB, SNR = 15dB, L=6
disp('alpha = 0.22, SNR = 15dB, SIR = 15dB, L = 6')
trunc = 5;

alpha = 0.22;
snr = 15;
SIR = 15;
L = 6;


BER_ISI_CCI_i_truncated('RC', alpha, snr, L, SIR, trunc)
BER_ISI_CCI_i_truncated('BTRC', alpha, snr, L, SIR, trunc)
BER_ISI_CCI_i_truncated('IPLCP', alpha, snr, L, SIR, trunc)
BER_ISI_CCI_i_truncated('ELP', alpha, snr, L, SIR, trunc)


disp('alpha = 0.35, SNR = 15dB, SIR = 15dB, L = 6')
alpha = 0.35;

BER_ISI_CCI_i_truncated('RC', alpha, snr, L, SIR, trunc)
BER_ISI_CCI_i_truncated('BTRC', alpha, snr, L, SIR, trunc)
BER_ISI_CCI_i_truncated('IPLCP', alpha, snr, L, SIR, trunc)
BER_ISI_CCI_i_truncated('ELP', alpha, snr, L, SIR, trunc)

disp('alpha = 0.5, SNR = 15dB, SIR = 15dB, L = 6')
alpha = 0.5;

BER_ISI_CCI_i_truncated('RC', alpha, snr, L, SIR, trunc)
BER_ISI_CCI_i_truncated('BTRC', alpha, snr, L, SIR, trunc)
BER_ISI_CCI_i_truncated('IPLCP', alpha, snr, L, SIR, trunc)
BER_ISI_CCI_i_truncated('ELP', alpha, snr, L, SIR, trunc)

trunc = 10;

disp('alpha = 0.22, SNR = 15dB, SIR = 15dB, L = 6, trunc=10')
BER_ISI_CCI_i_truncated('RC', alpha, snr, L, SIR, trunc)
BER_ISI_CCI_i_truncated('BTRC', alpha, snr, L, SIR, trunc)
BER_ISI_CCI_i_truncated('IPLCP', alpha, snr, L, SIR, trunc)
BER_ISI_CCI_i_truncated('ELP', alpha, snr, L, SIR, trunc)


disp('alpha = 0.35, SNR = 15dB, SIR = 15dB, L = 6')
alpha = 0.35;

BER_ISI_CCI_i_truncated('RC', alpha, snr, L, SIR, trunc)
BER_ISI_CCI_i_truncated('BTRC', alpha, snr, L, SIR, trunc)
BER_ISI_CCI_i_truncated('IPLCP', alpha, snr, L, SIR, trunc)
BER_ISI_CCI_i_truncated('ELP', alpha, snr, L, SIR, trunc)

disp('alpha = 0.5, SNR = 15dB, SIR = 15dB, L = 6')
alpha = 0.5;

BER_ISI_CCI_i_truncated('RC', alpha, snr, L, SIR, trunc)
BER_ISI_CCI_i_truncated('BTRC', alpha, snr, L, SIR, trunc)
BER_ISI_CCI_i_truncated('IPLCP', alpha, snr, L, SIR, trunc)
BER_ISI_CCI_i_truncated('ELP', alpha, snr, L, SIR, trunc)