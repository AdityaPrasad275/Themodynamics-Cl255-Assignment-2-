T = [-125 -130 -140 -150 -160 -170 -175 -180 -183]+273.15;
delH = zeros(1, size(T, 2));
delS = zeros(1, size(T, 2));
for i = 1:size(T, 2)
    delH(1, i) = delHS(T(1, i), 2);
    delS(1, i) = delHS(T(1, i), 1);
end
plot(T, delS);
xlabel('Temperature K');
ylabel('del S (Sv - Sl) (J/mol)');
figure;
plot(T, delH/1e3);
xlabel('Temperature K');
ylabel('del H (Hv - Hl) (kJ/mol)')
