%%execute all ISI

disp('alpha = 0.22, SNR = 10dB, trunc = 5')

alpha = 0.22;
snr = 10;
trunc = 5; %valor de truncado

BERi_truncated('RC', alpha, snr, trunc)
BERi_truncated('BTRC', alpha, snr, trunc)
BERi_truncated('IPLCP', alpha, snr, trunc)
BERi_truncated('ELP', alpha, snr, trunc)

disp('alpha = 0.35, SNR = 10dB, trunc = 5')

alpha = 0.35;
snr = 10;

BERi_truncated('RC', alpha, snr, trunc)
BERi_truncated('BTRC', alpha, snr, trunc)
BERi_truncated('IPLCP', alpha, snr, trunc)
BERi_truncated('ELP', alpha, snr, trunc)

alpha = 0.50;
snr = 10;



disp('alpha = 0.5, SNR = 10dB, trunc = 5')

BERi_truncated('RC', alpha, snr, trunc)
BERi_truncated('BTRC', alpha, snr, trunc)
BERi_truncated('IPLCP', alpha, snr, trunc)
BERi_truncated('ELP', alpha, snr, trunc)

disp('alpha = 0.22, SNR = 20dB, trunc = 5')

alpha = 0.22;
snr = 20;

BERi_truncated('RC', alpha, snr, trunc)
BERi_truncated('BTRC', alpha, snr, trunc)
BERi_truncated('IPLCP', alpha, snr, trunc)
BERi_truncated('ELP', alpha, snr, trunc)


disp('alpha = 0.35, SNR = 20dB, trunc = 5')


alpha = 0.35;
snr = 20;

BERi_truncated('RC', alpha, snr, trunc)
BERi_truncated('BTRC', alpha, snr, trunc)
BERi_truncated('IPLCP', alpha, snr, trunc)
BERi_truncated('ELP', alpha, snr, trunc)

alpha = 0.50;
snr = 20;

disp('alpha = 0.5, SNR = 20dB')

BERi_truncated('RC', alpha, snr, trunc)
BERi_truncated('BTRC', alpha, snr, trunc)
BERi_truncated('IPLCP', alpha, snr, trunc)
BERi_truncated('ELP', alpha, snr, trunc)

trunc = 10;
alpha = 0.22;
snr = 10;
disp('alpha = 0.22, SNR = 10dB')
disp('trunc = 10')

BERi_truncated('RC', alpha, snr, trunc)
BERi_truncated('BTRC', alpha, snr, trunc)
BERi_truncated('IPLCP', alpha, snr, trunc)
BERi_truncated('ELP', alpha, snr, trunc)

disp('alpha = 0.35, SNR = 10dB')

alpha = 0.35;
snr = 10;

BERi_truncated('RC', alpha, snr, trunc)
BERi_truncated('BTRC', alpha, snr, trunc)
BERi_truncated('IPLCP', alpha, snr, trunc)
BERi_truncated('ELP', alpha, snr, trunc)

alpha = 0.50;
snr = 10;



disp('alpha = 0.5, SNR = 10dB')

BERi_truncated('RC', alpha, snr, trunc)
BERi_truncated('BTRC', alpha, snr, trunc)
BERi_truncated('IPLCP', alpha, snr, trunc)
BERi_truncated('ELP', alpha, snr, trunc)

disp('alpha = 0.22, SNR = 20dB')

alpha = 0.22;
snr = 20;

BERi_truncated('RC', alpha, snr, trunc)
BERi_truncated('BTRC', alpha, snr, trunc)
BERi_truncated('IPLCP', alpha, snr, trunc)
BERi_truncated('ELP', alpha, snr, trunc)


disp('alpha = 0.35, SNR = 20dB')


alpha = 0.35;
snr = 20;

BERi_truncated('RC', alpha, snr, trunc)
BERi_truncated('BTRC', alpha, snr, trunc)
BERi_truncated('IPLCP', alpha, snr, trunc)
BERi_truncated('ELP', alpha, snr, trunc)

alpha = 0.50;
snr = 20;

disp('alpha = 0.5, SNR = 20dB')

BERi_truncated('RC', alpha, snr, trunc)
BERi_truncated('BTRC', alpha, snr, trunc)
BERi_truncated('IPLCP', alpha, snr, trunc)
BERi_truncated('ELP', alpha, snr, trunc)