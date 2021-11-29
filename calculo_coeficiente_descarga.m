close all
clear all
clc

%Datos constructivos
beta = 0.4;
D = 3*25.4*1e-3; %[m]
d = D*beta; %[m]

%Datos medicion 1
dP1 = [190 240 390 610 714 850 980 1100 1220 1350 1450 1570 1800 1900 2000 2120 2220 2340 2465]; %presion diferencial [Pa] medida en la garganta
Qmed1 = 1.0.*[36.8 53 64.3 81.5 91 99 105.5 111.5 117 122.3 125.5 130.6 144.5 148.3 152.5 156.5 160.2 164 168]; %caudal volumetrico medido con el flujometro[m^3/h]
Qmed1 = Qmed1./3600; %[m3/s]
T1 = [25.5 26.5 27 27.6 28.3 28.5 28.9 29.3 30 30.4 30.8 30.9 31.2 31.5 31.5 31.7 31.8 31.9 32.1]; %temperatura [°C]
Hr1 = [57 52.2 51.2 50.2 48.9 48.1 47.3 45.9 44.1 42.9 42.3 42.3 41.1 41.4 40.7 39.9 38.4 39.3 39.2]; %humedad relativa [%]
Vp1 = [3.61 3.618 3.617 3.615 3.617 3.615 3.61 3.59 3.59 3.6 3.62 3.61 3.60 3.6 3.59 3.58 3.59 3.6 3.59]; %tension para presion absoluta [V]. R = 469 ohm.

%Datos medicion 2
dP2 = [2450 2350 2220 2130 2020 1900 1890 1680 1570 1460 1338 1225 1100 980 850 710 615 385 248 224];
Qmed2 = 1.0.*[167.5 164 160 156 152 148.5 144.6 140.3 130.3 126.3 121.3 116.3 111.3 105.3 98.5 90.5 81.2 63.5 53 36.8];
Qmed2 = Qmed2./3600;
T2 = [30.7 31.4 31.6 31.9 32.2 32.3 32.5 32.5 32.7 33.3 33.6 33.8 34 34 33.5 33.4 33 33.3 33.3 33.1];
Hr2 = [44.4 40.2 39.5 39.1 39.2 38.6 38.2 38 37.7 36.8 36 35.7 35.7 35.7 36.1 37.1 37.1 36.6 36.8 37];
Vp2 = [3.6 3.61 3.61 3.6 3.59 3.58 3.57 3.56 3.54 3.56 3.56 3.57 3.59 3.59 3.59 3.6 3.6 3.6 3.6 3.6 ];

%Datos totales
Qmed = [Qmed1 Qmed2];
[Qmed,ord] = sort(Qmed);

dP = [dP1 dP2];
dP = dP(ord);
T = [T1 T2];
T = T(ord);
Hr = [Hr1 Hr2];
Hr = Hr(ord);
Vp = [Vp1 Vp2];
Vp = Vp(ord);


Ro = zeros(1,length(T));
Roa = zeros(1,length(T));
Mu = zeros(1,length(T));
I = zeros(1,length(T));
pressure = zeros(1,length(T));
eps = zeros(1,length(T));


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

n = 300;
Cd = linspace(0.93,0.95,n);
Qteo = zeros(n,length(Qmed)); %Una fila de Qteo para cada Cd

for i = 1:n
    
    Qteo(i,:) = Cd(i) * eps * (pi/4) * d^2 .* sqrt(2.*dP./(Roa.*(1-beta^4)));%[m3/s] %Caudal volumétrico calculado con los dP medidos
        
end

% plot(dP,Qmed,'bo')

clear i

SL1 = zeros(1,n);
sumSL1 = zeros(1,length(Qmed));
SL2 = zeros(1,n);
sumSL2 = zeros(1,length(Qmed));

for j = 1:n
    
    for i = 1:length(Qmed)
        
        sumSL1(i) = abs(Qmed(i)-Qteo(j,i));
        sumSL2(i) = (Qmed(i)-Qteo(j,i))^2;

    end
    
    SL1(j) = sum(sumSL1);
    SL2(j) = sum(sumSL2);
end

%
figure
yyaxis left
plot(Cd,SL1,'ro')
ylabel('L1-norm')

hold on

yyaxis right
plot(Cd,SL2,'go')
xlabel('Coeficiente de descarga')
ylabel('L2-norm')

legend('SL1','SL2','Location','North')
title('L-norm vs Cd')
grid minor
%}

Cdfinal1 = Cd(find(SL1 == min(SL1)))
Cdfinal2 = Cd(find(SL2 == min(SL2)))

Qteo1 = Cdfinal1 * eps * (pi/4) * d^2 .* sqrt(2.*dP./(Roa.*(1-beta^4)));%[m3/s];
Qteo2 = Cdfinal2 * eps * (pi/4) * d^2 .* sqrt(2.*dP./(Roa.*(1-beta^4)));%[m3/s];

err1 = abs(Qteo1 - Qmed);
err2 = abs(Qteo2 - Qmed);


%Caudal másico
Qmteo1 = Cdfinal1 * eps * (pi/4) * d^2 .* sqrt(2.*Roa.*dP./((1-beta^4)));%[kg/s];;
Qmteo2 = Cdfinal2 * eps * (pi/4) * d^2 .* sqrt(2.*Roa.*dP./((1-beta^4)));%[kg/s];;

QmMed = Roa.*Qmed;

%plot Qvolumetrico
figure(2)

plot(dP,Qmed,'bo-')%,dP,Qteo1,'ro-',dP,Qteo2,'go-')
xlabel('Presión diferencial [Pa]')
ylabel('Caudal volumétrico [m^3/s]')
hold on
% errorbar(dP,Qteo1,err1,'ro-')
errorbar(dP,Qteo2,err2,'go-')
% legend('Qmedido','Qteorico1 (L1-norm)','Qteorico2 (L2-norm)','Location','Best')
legend('Q_{medido}','Q_{teórico}','Location','Best')
title('Caudal volumétrico vs presión diferencial')
grid minor

%plot Qmasico
figure

plot(dP,QmMed,'bo-',dP,Qmteo2,'go-')%,dP,Qmteo1,'ro-')
xlabel('Presión diferencial [Pa]')
ylabel('Caudal másico [kg/s]')
hold on
% % errorbar(dP,Qteo1,err1,'ro-')
% errorbar(dP,Qteo2,err2,'go-')
% legend('Qmedido','Qteorico1 (L1-norm)','Qteorico2 (L2-norm)','Location','Best')
legend('Qm_{medido}','Qm_{teórico}','Location','Best')
title('Caudal másico vs presión diferencial')
grid minor



%
%
n1=3;
ep1 = (err1./Qmed)*100;
ep2 = (err2./Qmed)*100;
mEp1 = mean(ep1(n1:end))*ones(1,length(dP));
mEp2 = mean(ep2(n1:end))*ones(1,length(dP));

figure
plot(dP(n1:end),ep1(n1:end),'ro-',dP(n1:end),ep2(n1:end),'go-')
hold on
plot(dP(n1:end),mEp1(n1:end),'r--',dP(n1:end),mEp2(n1:end),'g--')%,'LineWidth',2)
ylabel('\DeltaQ_v [%]')
xlabel('Presion diferencial [Pa]')
title('Diferencia porcentual entre Q_{teorico} y Q_{medido} vs Presion diferencial')
legend(['L1-norm, C_d = ',num2str(Cdfinal1)],['L2-norm, C_d = ',num2str(Cdfinal2)],...
    ['\DeltaQ_v promedio = ',num2str(mEp1(1)),'%'],['\DeltaQ_v promedio = ',num2str(mEp2(1)),'%'])
grid minor

%}
%Resultados de Cd sin la correccion del 5% del error de medicion
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
%}
