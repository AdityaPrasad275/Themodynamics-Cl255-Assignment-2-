function y = f(Z, T, P)

[A B] = pt_consts(T, P);

term1 = Z-1;
term2 = log(Z-B);
term3 = (Z + (1+sqrt(2))*B)/(Z + (1-sqrt(2))*B);
term3 = (A/(sqrt(8)*B))*log(term3);

y = exp(term1 - term2 - term3).*P;
end