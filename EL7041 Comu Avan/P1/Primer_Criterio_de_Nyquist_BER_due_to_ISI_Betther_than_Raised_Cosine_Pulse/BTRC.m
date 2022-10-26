function[v]=BTRC(t,alpha)
% Function for the BTRC pulse shape.

% Set the parameters
T=1;
beta = (pi * alpha) / log(2);

% Calculate the sinc expression
sinc_BTRC = sinc(t / T);

% Calculate the second expression for BTRC
numerador_second_BTRC = (2 * beta * (t / T)) * sin(pi * alpha * (t / T));
numerador_second_BTRC = numerador_second_BTRC + 2 * cos(pi * alpha * (t / T)) - 1;

% Never 0
denominador_second_BTRC = (beta * (t / T))^2 + 1;

% Get the BTRC Pulse
v = sinc_BTRC * numerador_second_BTRC / denominador_second_BTRC;

end