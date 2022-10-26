function[v]=LCP(t,alpha)
% Function for the BTRC pulse shape.

% Set the parameters
T=1;

%LINEAR COMBINATION PULSE
%alpha*PLP(n=1) +(1-alpha)RC
constante=1; %alpha

Exp1_Num=((1-constante)*sin(pi*(t / T))*cos(pi*alpha*(t / T))); %RC
Exp1_Den=pi*(t / T)*(1-(2*alpha*(t / T))^2);
%Exp1DenZero=find(abs(Exp1_Den) < 10^-10);
Exp1_Total=Exp1_Num./Exp1_Den;
%Exp1_Total(Exp1DenZero)=0.5;

Exp2_Num=constante*sin(pi*(t / T))*sin(pi*alpha*(t / T)); %PLP
Exp2_Den=pi^2*alpha*(t / T)^2;
%Exp2DenZero=find(abs(Exp2_Den) < 10^-10);
Exp2_Total=Exp2_Num./Exp2_Den;
%Exp2_Total(Exp2DenZero)=0.5;

Linear_Combination_Pulse_PLP=(Exp1_Total+Exp2_Total);

% Get the BTRC Pulse
v = Linear_Combination_Pulse_PLP;

end



