function[v]=BTRC_Ten(t,alpha)
% Function for the BTRC pulse shape truncated to 10 periods.

T=1;
% Truncate the values.
if t > 10
    v = 0;
else
    v = BTRC(t, alpha);
end
end