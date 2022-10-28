R = 8.314;
Tc=154.6; %K
Pc=5.046e6; %Pa
Vc=0.0732e-3; %m^3/mol
a=0.45724*(R*Tc)^2/Pc; %at T=Tc
b=0.0778*R*Tc/Pc;
omega=0.021;

T = [-183 -180 -175 -170 -160 -150 -140 -130 -100 -75 -50 -25 0 25 50 75 100 125] + 273.15;
P=logspace(0, log(80)/log(10),300)*1e5; %Pa

nT = size(T, 2);
nP = size(P, 2);

enth = zeros(1, nP);
for t = 1:nT
    for p = 1:nP
        enth(1, p) = 0.001*enthalpy(T(1, t), P(1, p));
    end
    plot(enth, P/1e5);
    hold on;
end
xlabel('Enthalpy (kJ/mol)');
ylabel('Pressure (bar)');
figure;

T = 73.15:10:400.15;
P = [1 2 5 10 20 50]*1e5;
nT = size(T, 2);
nP = size(P, 2);
ent = zeros(1, nT);
for p1 = 1:nP
    for t1 = 1:nT
        ent(1, t1) = entropy(T(1, t1), P(1, p1));
    end
    plot(ent, T);
    hold on;
end
xlabel('Entropy (J/mol)');
ylabel('Temperature (K)');
legend('1 bar', '2 bar', '5 bar', '10 bar', '20 bar', '50 bar');
figure

T = [-183 -180 -175 -170 -160 -150 -140 -130 -100 -75 -50 -25 0 25 50 75 100 125] + 273.15;
P=logspace(0,log(80)/log(10),300)*1e5; %Pa
nT = size(T, 2);
nP = size(P, 2);
ent = zeros(1, nP);
for t1 = 1:nT
    for p1 = 1:nP
        ent(1, p1) = entropy(T(1, t1), P(1, p1));
    end
    plot(ent, P/1e5);
    hold on;
end
xlabel('Entropy (J/mol)');
ylabel('Pressure (bar)');