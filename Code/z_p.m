Pc = 5.046e6;
Tc = 154.6;
R = 8.314;
b = 0.0778*(R*Tc)/Pc;
a = 0.45724*((R*Tc)^2)/Pc;

Pr = logspace(log(0.1)/log(10), log(7)/log(10), 300);
P = Pr*Pc;

Tr = [1 1.1 1.2 1.3 1.5 2];
T = Tr*Tc;

z = zeros(1, size(P, 2));

for i = 1:size(T, 2)
    for j = 1:size(P, 2)
        [A, B] = pt_consts(T(1, i), P(1, j));
        alpha= -1.+B;
        beta=A-3*B*B-2*B;
        gamma= -A*B + B*B*(B+1.);

        PReq=[1.0 alpha beta gamma];
        zr=real(roots(PReq));

        %get only relevant values of z
        check_zr=abs(polyval(PReq,zr))<1.e-8;
        pos=find(check_zr==1);
        z(1, j) = zr(pos);
    end
    plot(Pr, z);
    xlabel('Reduced Pressure');
    ylabel('Compressiblity Factor Z  (PV/RT)');
    legend('Tr = 1', 'Tr = 1.1', 'Tr = 1.2', 'Tr = 1.3', 'Tr = 1.5', 'Tr = 2');
    hold on
end
