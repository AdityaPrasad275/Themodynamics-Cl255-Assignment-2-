function z=PR(T,P)

  [A, B] = pt_consts(T, P);
  alpha= -1.+B;
  beta=A-3*B*B-2*B;
  gamma= -A*B + B*B*(B+1.);

  PReq=[1.0 alpha beta gamma];
  zr=real(roots(PReq));

  %get only relevant values of z
  check_zr=abs(polyval(PReq,zr))<1.e-8;
  pos=find(check_zr==1);
  z=zr(pos);
end