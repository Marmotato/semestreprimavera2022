function[v]=IPLCP(t,alpha)

% IMPROVED LINEAR COMBINATION PULSE
% ((mu*A +(1-mu)B) * sinc * exp)^gamma

T=1;

mu = 1.60;
eps = 0.1;
gamma = 1;

sincNum = sin(pi*(t / T)); % numerator of the sinc function
sincDen = (pi*(t / T)); % denominator of the sinc function
sincOp = sincNum./sincDen;


ExpIPLC = exp(-eps*pi^2*((t / T)).^2);

IPLCNum = 4*(1 - mu)*(sin(alpha*pi*(t / T)/2).^2) + alpha*pi*mu*(t / T).*sin(alpha*pi*(t / T));
IPLCDen = (alpha*pi*(t / T)).^2;
IPLC_1 = IPLCNum ./ IPLCDen;


% Get the BTRC Pulse
v = ExpIPLC.*(sincOp.*IPLC_1).^gamma;



end

