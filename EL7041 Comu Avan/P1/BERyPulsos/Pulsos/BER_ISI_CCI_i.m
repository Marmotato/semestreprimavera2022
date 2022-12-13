function [BER]=BER_ISI_CCI_i(str,alpha,snr, L, SIRdB)
%Alpha: Factor de roll-off
%L: Elementos generadores de CCI
%STRING --> NOMBRE DEL PULSO --> CARPETA CON PULSOS
%BER: Valores del Bit-Error rate para cada offset
%Cada columna corresponde a los offsets t/T [0.05, 0.1, 0.2, 0.25]
%
%% Parametros de Simulacion 

% Set the parameters for the simulation
T=1;
nbits = 2^10;
N     = floor(nbits/2);     
M     = 100;                
omega = 0.10; 



% Obtain the value for the SNR and SIR in linear
coeff = 10^(snr/20);    
coeff2 = 10^(SIRdB/20);

ri = sqrt(1/L) * (coeff/coeff2);               
  

% Calculate the offsets
offset = [0.05, 0.1, 0.2 0.25];

% Limits to use for the sum
a=linspace(-N,-1,N);
b=linspace(1,N,N);
ab=[a b];
values = [-1, 1];

% Pulse shape to use
fh=str2func(str);

% Calculate the BER
%cd('Pulsos')
sumaT2=zeros(1,length(offset));
for c=1:length(offset)
    g0 = coeff * fh(offset(c) * T,alpha);
    gk=zeros(length(ab),1);
    
    % Calculate the values for gk
    for i=1:length(ab)
    end
    
    % Calculate the sum and product
    suma=0;
    mult1=1;
    mult2=1;
    % Calculate the sum
    for m=1:2:M
        % Calculate the products
        %CCI
        for i=1:L
            mult1 = double(mult1 *besselj(0, (m*omega*ri ) ) ); 
        end
        %ISI
        for k=1:length(gk)
            mult2= double(mult2 * cos(m*omega*gk(k)));
        end
        suma= double(suma + ((exp(-(m * omega)^2 / 2) * sin(m * omega * g0))/m) * (mult1 * mult2));
        mult1=1;
        mult2=1;
    end
    sumaT2(c) = (1./2. - (2./pi) * suma);
end

disp(str)

BER=sumaT2;
%BER=vpa(BER,8); %Precisión del BER

end