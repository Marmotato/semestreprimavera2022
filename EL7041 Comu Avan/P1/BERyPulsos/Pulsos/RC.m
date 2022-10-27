function[v]=RC(t,alpha)
% Function for the BTRC pulse shape.

% Set the parameters
T=1;


sincNum = sin(pi*(t / T)); % numerator of the sinc function
sincDen = (pi*(t / T)); % denominator of the sinc function
sincOp = sincNum./sincDen;


cosNum = cos(alpha*pi*(t / T));
cosDen = (1-(2*alpha*(t / T)).^2);
cosOp = cosNum./cosDen;


RC_Filter = sincOp.*cosOp;

% Get the BTRC Pulse
v = RC_Filter;

end