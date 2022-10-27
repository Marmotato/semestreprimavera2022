%%execute all CCI and ISI

%% SIR = 15dB, SNR = 15dB, L=6
disp('alpha = 0.22, SNR = 15dB, SIR = 15dB, L = 6')

alpha = 0.22;
snr = 15;
SIR = 15;
L = 6;


BER_ISI_CCI_i('RC', alpha, snr, SIR, L)
BER_ISI_CCI_i('BTRC', alpha, snr, SIR, L)
BER_ISI_CCI_i('IPLCP', alpha, snr, SIR, L)
BER_ISI_CCI_i('ELP', alpha, snr, SIR, L)


disp('alpha = 0.35, SNR = 15dB, SIR = 15dB, L = 6')
alpha = 0.35;

BER_ISI_CCI_i('RC', alpha, snr, SIR, L)
BER_ISI_CCI_i('BTRC', alpha, snr, SIR, L)
BER_ISI_CCI_i('IPLCP', alpha, snr, SIR, L)
BER_ISI_CCI_i('ELP', alpha, snr, SIR, L)

disp('alpha = 0.5, SNR = 15dB, SIR = 15dB, L = 6')
alpha = 0.5;

BER_ISI_CCI_i('RC', alpha, snr, SIR, L)
BER_ISI_CCI_i('BTRC', alpha, snr, SIR, L)
BER_ISI_CCI_i('IPLCP', alpha, snr, SIR, L)
BER_ISI_CCI_i('ELP', alpha, snr, SIR, L)
