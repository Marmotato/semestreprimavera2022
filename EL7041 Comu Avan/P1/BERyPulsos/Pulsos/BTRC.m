function[v]=BTRC(t,alpha)
% Function for the BTRC pulse shape.

% Set the parameters
T=1;
beta2 = (pi * alpha) / log(2);

% Calculate the sinc expression
sinc_BTRC_num = sin(pi.*(t / T));
sinc_BTRC_den = pi.*(t / T);
sinc_BTRC_tot = sinc_BTRC_num./sinc_BTRC_den;

% Calculate the second expression for BTRC
numerador_second_BTRC = (2 * beta2 * (t / T)) .* sin(pi * alpha * (t / T));
numerador_second_BTRC = numerador_second_BTRC + (2 .* cos(pi * alpha * (t / T)) - 1);

% Never 0
denominador_second_BTRC = (beta2 * (t / T)).^2 + 1;

% Get the BTRC Pulse
v = ((sinc_BTRC_tot .* numerador_second_BTRC) ./ denominador_second_BTRC);
v(isnan(v))=1; %corregimos en el punto 0

end