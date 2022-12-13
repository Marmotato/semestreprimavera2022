function[v]=Truncate(t,alpha, str, trunc)
% Function for the BTRC pulse shape truncated to "trunc" periods.
T=1;

fh=str2func(str);

tabs = abs(t);

% Truncate the value
if tabs > trunc
    v = 0;
else
    v = fh(t, alpha);
end
end