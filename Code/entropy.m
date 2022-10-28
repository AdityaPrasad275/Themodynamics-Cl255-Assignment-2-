function y = entropy(T, P)
[A, B] = pt_consts(T, P);
alpha= -1.+B;
beta=A-3*B*B-2*B;
gamma= -A*B + B*B*(B+1.);
PReq=[1.0 alpha beta gamma];
zr=real(roots(PReq));
if T < 173.15
    if P< p_vap(T)
        Z = max(zr);
    elseif P>p_vap(T)
        Z = min(zr);
    else
        return
    end
else 
    Z = PR(T, P);
end

Pc = 5.046e6;
Tc = 154.6;
R = 8.314;
b = 0.0778*(R*Tc)/Pc;
a = 0.45724*((R*Tc)^2)/Pc;

kappa = 0.4069;
alpha=1.+kappa*(1.-sqrt(T/Tc));
da_dt = (-1)*a*alpha*(kappa/sqrt(T*Tc));
c1 = (da_dt)/(sqrt(8)*b);
B = (b*P)/(R*T);
c2 = (Z + (1+sqrt(2))*B)/(Z + (1-sqrt(2))*B);
c2 = log(c2);
s1 = R*log(Z-B) + c1*c2;
c3 = 25.46*log(T/298.15);
c4 = (0.01519)*(T - 298.15);
c5 = ((0.7151e-5)/2)*(T^2 - 298.15^2);
c6 = ((1.311e-9)/3)*(T^3-298.15^3);
s2 = c3 + c4 - c5 + c6 - R*log(P/1e5);
y = s1 + s2;
end