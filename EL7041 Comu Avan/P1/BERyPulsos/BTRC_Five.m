function[v]=BTRC_Five(t,alpha)
% Function for the BTRC pulse shape truncated to 5 periods.
T=1;

% Truncate the value
if t > 5
    v = 0;
else
    v = BTRC(t, alpha);
end
end