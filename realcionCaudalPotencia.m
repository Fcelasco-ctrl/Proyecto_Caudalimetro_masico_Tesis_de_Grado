clear
home

t = 20; %[°C]
P = 101325 + 3000;%[Pa]
H = 20; %[%]
[ro,roa] = density(t,H,P);

LHV = 43.4; %[MJ/kg]
A_C = 14.7/1; %Kg_aire/Kg_gasolina
Qv_max = 0.11; %[m^3/s]
Eta = 0.37; %Eficiencia
Pot = 150; %[kW] Potencia del motor

Qm_max = roa*Qv_max; %Caudal masico de aire maximo [kg/s]

Qm_comb = Qm_max/A_C; %Caudal masico de combustible [kg/s]

Pmax = Eta*Qm_comb*LHV*1000; %[kW]




