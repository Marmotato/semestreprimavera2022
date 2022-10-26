function[v]=RC(t,alpha)
% Function for the BTRC pulse shape.

% Set the parameters
T=1;

sincNum = sin(pi*(t / T)); % numerator of the sinc function
sincDen = (pi*(t / T)); % denominator of the sinc function
%sincDenZero = find(abs(sincDen) < 10^-10);
sincOp = sincNum/sincDen;
%sincOp(sincDenZero) = 1; % sin(pix/(pix) =1 for x =0

cosNum = cos(alpha*pi*(t / T));
cosDen = (1-(2*alpha*(t / T))^2);
%cosDenZero = find(abs(cosDen)<10^-10);
cosOp = cosNum/cosDen;
%cosOp(cosDenZero) = pi/4;

RC_Filter = sincOp*cosOp;

% Get the BTRC Pulse
v = RC_Filter;

end