function [A_z, B_z] = pt_consts(T, P)
Pc = 5.046e6;
Tc = 154.6;
R = 8.314;
b = 0.0778*(R*Tc)/Pc;
a = 0.45724*((R*Tc)^2)/Pc;
 
kappa = 0.4069;
alpha=1+kappa*(1-sqrt(T/Tc));
alpha=alpha*alpha;
a=a*alpha;

A_z = (a*P)/((R*T)^2);
B_z = (b*P)/(R*T);

end