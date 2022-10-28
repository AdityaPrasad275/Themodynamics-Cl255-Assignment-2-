R = 8.314;
v = logspace(-2, 1, 500)*1e-3;
p = zeros(1, size(v, 2));
t = [-175 -150 -125]+273.15;

for j =1:size(t, 2)
    vapour_p = p_vap(t(1, j));
    [A, B] = pt_consts(t(1, j), vapour_p);
    alpha= -1.+B;
    beta=A-3*B*B-2*B;
    gamma= -A*B + B*B*(B+1.);

    PReq=[1.0 alpha beta gamma];
    zr=real(roots(PReq));
    check_zr=abs(polyval(PReq,zr))<1.e-8 & zr>1.e-5;
   
    zl = zr(3);
    zv = zr(1);

    vl = zl*R*t(1, j)/vapour_p;
    vu = zv*R*t(1, j)/vapour_p;

    for i = 1:size(v, 2)
        if v(1, i) <vl
            p(1, i) = pv(v(1, i), t(1, j));
        elseif v(1, i) > vu
            p(1, i) = pv(v(1, i), t(1, j));
        else
            p(1, i) = vapour_p;
        end
    end
    
    semilogx(v, p);
    xlim([2e-5 0.01]);
    ylim([1e5 100e5]);
    hold on;
end

xlabel('Volume (m^3/mol)');
ylabel('Pressure (Pa)');
legend('T = -175 C', 'T = -150 C', 'T = -125 C');
