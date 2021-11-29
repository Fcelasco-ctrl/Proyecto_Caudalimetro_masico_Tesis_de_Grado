close all
clear all
clc

%Datos constructivos
beta = 0.4;
D = 3*25.4*1e-3; %[m]
d = D*beta; %[m]

%Datos medicion 1
dP = [9 12 26 13 18 45]; %presion diferencial [Pa] medida en la garganta
T = [26.5 26.4 26.4 26.3 26.4 26.4]; %temperatura [°C]
Hr = [54 54.2 54.5 54.7 54.5 54.9]; %humedad relativa [%]
Vp = [3.714 3.68 3.7 3.69 3.7 3.69]; %tension para presion absoluta [V]. R = 469 ohm.
rpm =  [1000 1000 1500 1500 2000 2000]; %RPM del motor
tita = [35 80 80 35 35 80]; %angulo de mariposa [°]

Ro = zeros(1,length(T));
Roa = zeros(1,length(T));
Mu = zeros(1,length(T));
I = zeros(1,length(T));
pressure = zeros(1,length(T));
eps = zeros(1,length(T));
% Cd = 0.9408;
Cd = 0.9343;
for i = 1:length(T)
    %Presión
    
    R = 469; %resistencia [ohm]
    pmax = 2.5; %pmax del rango sensor de presion [bar]
    imin = 4e-3; %corriente minima de salida del sensor [A]
    imax = 20e-3; %corriente maxima de salida del sensor [A]
    
    I(i) = Vp(i)/R; %corriente equivalente de presion
    
    pressure(i) = (I(i)/imax)*pmax*1e5; %[Pa]
    
    
    
    %Densidad
    
    Pamb = pressure(i);
    P = 0;
    t = T(i);
    Hum = Hr(i);
    Tamb = 26;%[°C]
    
    [ro,roa] = density(t,Hum,P,Pamb); %ro:densidad aire humedo -  roa:densidad aire seco
    [rho,mu] = moistAir2(Tamb,Hum,P,Pamb,t);%rho:no se usa - mu:viscosidad
 
    Ro(i) = ro; %densidad aire húmedo [kg/m3]
    Roa(i) = roa; %densidad aire seco [kg/m3]
    Mu(i) = mu; %viscosidad
    
    
    P2 = pressure(i) - dP(i);
    k = 1.4;
    r = P2./pressure(i);
    eps(i) =  sqrt((r^(2/k))*(k/(k-1))*((1-r^((k-1)/k))/(1-r))*((1-beta^4)/(1-(beta^4)*(r^(2/k)))));% [dimlss] Factor de correccion por expansion. (Fluidos compresibles)


   
end

clear i

Qteo = zeros(1,length(T)); %Una fila de Qteo

for i = 1:length(T)
    
    Qteo(i) = Cd * eps(i) * (pi/4) * d^2 .* sqrt((2.*dP(i)*Roa(i))/((1-beta^4)));%[m3/s] %Caudal másico calculado con los dP medidos
        
end


ord35 = find(tita == 35);
ord80 = find(tita == 80);


figure
plot(rpm(ord35),Qteo(ord35),'bo-',rpm(ord80),Qteo(ord80),'rs-')
xlabel('Velocidad de motor [r.p.m.]')
ylabel('Caudal másico [kg/s]')
legend('\theta = 35°','\theta = 80°','Location','Best')
title('Caudal másico vs Velocidad de motor')
grid minor




%{
ep1 = (err1./Qmed)*100;
ep2 = (err2./Qmed)*100;
figure
plot(dP(3:end),ep1(3:end),'ro-',dP(3:end),ep2(3:end),'go-')
ylabel('Delta Qv [%]')
xlabel('Presion diferencial [Pa]')
title('Diferencia porcentual entre Qteo y Qmed vs Presion diferencial')
legend('Qteorico1 (L1-norm)','Qteorico2 (L2-norm)')
grid minor

%}

%Datos totales
% Cdfinal1 = 0.9408
% Cdfinal2 = 0.9343

%Datos nro 1
% Cdfinal1 = 0.9456
% Cdfinal2 =  0.9386

%Datos nro 2
% Cdfinal1 = 0.9315
% Cdfinal2 =0.9304

%{
dPc = 1e4*[    0.0108    0.0134    0.0162    0.0194    0.0228    0.0265    0.0305    0.0348    0.0394    0.0442    0.0493    0.0547...    0.0604,
    0.0663    0.0725    0.0791    0.0858    0.0929    0.1003    0.1079    0.1158    0.1240    0.1324    0.1412    0.1502    0.1595...
    0.1691    0.1790    0.1891    0.1995    0.2102    0.2212    0.2325    0.2440    0.2559    0.2680    0.2803    0.2930    0.3059...
    0.3192    0.3327    0.3464    0.3605    0.3748    0.3895    0.4044    0.4195    0.4350    0.4507    0.4667    0.4830    0.4996...
    0.5165    0.5336    0.5510    0.5687    0.5867    0.6049    0.6235    0.6423    0.6614    0.6808    0.7004    0.7203    0.7406...
    0.7610    0.7818    0.8029    0.8242    0.8458    0.8677    0.8899    0.9123    0.9350    0.9580    0.9813    1.0049    1.0288...
    1.0529    1.0773];
Qteorico1 = Cdfinal1 * mean(eps) * (pi/4) * d^2 .* sqrt(2.*dPc./(mean(Roa).*(1-beta^4)));%[m3/s];
Qteorico2 = Cdfinal2 * mean(eps) * (pi/4) * d^2 .* sqrt(2.*dPc./(mean(Roa).*(1-beta^4)));%[m3/s];

figure
% plot(dP,Qmed,'bo-')%,dP,Qteo1,'ro-',dP,Qteo2,'go-')
xlabel('Presion diferencial [Pa]')
ylabel('Caudal volumétrico [m^3/s]')
hold on
errorbar(dP,Qteo1,err1,'ro-')
errorbar(dP,Qteo2,err2,'go-')
plot(dPc,Qteorico1,dPc,Qteorico2)
legend('Qmedido','Qteorico1 (L1-norm)','Qteorico2 (L2-norm)','Location','Best')
title('Caudal volumétrico vs presion diferencial')
grid minor
axis([0 2500 0.0099 0.05])
%}