%%execute all ISI

disp('alpha = 0.22, SNR = 10dB')

alpha = 0.22;
snr = 10;

BERi('RC', alpha, snr)
BERi('BTRC', alpha, snr)
BERi('IPLCP', alpha, snr)
BERi('ELP', alpha, snr)

disp('alpha = 0.35, SNR = 10dB')

alpha = 0.35;
snr = 10;

BERi('RC', alpha, snr)
BERi('BTRC', alpha, snr)
BERi('IPLCP', alpha, snr)
BERi('ELP', alpha, snr)

alpha = 0.50;
snr = 10;



disp('alpha = 0.5, SNR = 10dB')

BERi('RC', alpha, snr)
BERi('BTRC', alpha, snr)
BERi('IPLCP', alpha, snr)
BERi('ELP', alpha, snr)

disp('alpha = 0.22, SNR = 20dB')

alpha = 0.22;
snr = 20;

BERi('RC', alpha, snr)
BERi('BTRC', alpha, snr)
BERi('IPLCP', alpha, snr)
BERi('ELP', alpha, snr)


disp('alpha = 0.35, SNR = 20dB')


alpha = 0.35;
snr = 20;

BERi('RC', alpha, snr)
BERi('BTRC', alpha, snr)
BERi('IPLCP', alpha, snr)
BERi('ELP', alpha, snr)

alpha = 0.50;
snr = 20;

disp('alpha = 0.5, SNR = 20dB')

BERi('RC', alpha, snr)
BERi('BTRC', alpha, snr)
BERi('IPLCP', alpha, snr)
BERi('ELP', alpha, snr)