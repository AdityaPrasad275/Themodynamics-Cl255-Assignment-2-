function y = p_vap(T)
Pc = 5.046e6;
Pvap = Pc; %Pa -- initial guess
nroots=1;
while nroots~=3 %check whether this initial guess provides 3 roots
    z=PR(T,Pvap);
    nroots=size(z,1);
    if nroots==1
        Pvap=Pvap-10e5; %New guess
    end
end
fl = 1;
fv = 0.1;
while abs((fl/fv)-1) > 0.0001
    [A, B] = pt_consts(T, Pvap);
    alpha= -1.+B;
    beta=A-3*B*B-2*B;
    gamma= -A*B + B*B*(B+1.);

    PReq=[1.0 alpha beta gamma];
    zr=real(roots(PReq));
    check_zr=abs(polyval(PReq,zr))<1.e-8 & zr>1.e-5;
   
    zl = zr(3);
    zv = zr(1);

    fl = f(zl, T, Pvap);
    fv = f(zv, T, Pvap);
    Pvap = Pvap*fl/fv;
end
y = Pvap;
end