function[v]=ELP(t,alpha)
% Function for the BTRC pulse shape.

% Set the parameters
T=1;

%EXPONENTIAL LINEAR PULSE
constante=1; %alpha
beta = 0.1; %Beta


Exp2_2_Num=constante*sin(pi*(t / T))*sin(pi*alpha*(t / T)); %PLP
Exp2_2_Den=pi^2*alpha*(t / T)^2;
Exp2_2_Total=Exp2_2_Num/Exp2_2_Den;

Exp3 = exp(-pi*(beta/2)*(t / T)^2);

Exponential_Linear_Pulse=(Exp3.*Exp2_2_Total);

% Get the BTRC Pulse
v = Exponential_Linear_Pulse;

end


