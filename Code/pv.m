function y = pv(V, T)
Pc = 5.046e6;
Tc = 154.6;
R = 8.314;
b = 0.0778*(R*Tc)/Pc;
a = 0.45724*((R*Tc)^2)/Pc;

kappa = 0.4069;
alpha=1+kappa*(1-sqrt(T/Tc));
alpha=alpha*alpha;
a=a*alpha;
term1 = (R*T)/(V-b);
term2 = a/(V^2 + 2*V*b - b^2);
y = term1 - term2;
end