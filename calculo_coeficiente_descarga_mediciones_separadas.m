close all
clear all
clc

%Datos constructivos
beta = 0.4;
D = 3*25.4*1e-3; %[m]
d = D*beta; %[m]

%Datos medicion 1
dP1 = [190 240 390 610 714 850 980 1100 1220 1350 1450 1570 1800 1900 2000 2120 2220 2340 2465]; %presion diferencial [Pa] medida en la garganta
Qmed1 = 1.05*[36.8 53 64.3 81.5 91 99 105.5 111.5 117 122.3 125.5 130.6 144.5 148.3 152.5 156.5 160.2 164 168]; %caudal volumetrico medido con el flujometro[m^3/h]
Qmed1 = Qmed1./3600; %[m3/s]
T1 = [25.5 26.5 27 27.6 28.3 28.5 28.9 29.3 30 30.4 30.8 30.9 31.2 31.5 31.5 31.7 31.8 31.9 32.1]; %temperatura [°C]
Hr1 = [57 52.2 51.2 50.2 48.9 48.1 47.3 45.9 44.1 42.9 42.3 42.3 41.1 41.4 40.7 39.9 38.4 39.3 39.2]; %humedad relativa [%]
Vp1 = [3.61 3.618 3.617 3.615 3.617 3.615 3.61 3.59 3.59 3.6 3.62 3.61 3.60 3.6 3.59 3.58 3.59 3.6 3.59]; %tension para presion absoluta [V]. R = 469 ohm.

%Datos medicion 2
dP2 = [2450 2350 2220 2130 2020 1900 1890 1680 1570 1460 1338 1225 1100 980 850 710 615 385 248 224];
Qmed2 = 1.05*[167.5 164 160 156 152 148.5 144.6 140.3 130.3 126.3 121.3 116.3 111.3 105.3 98.5 90.5 81.2 63.5 53 36.8];
Qmed2 = Qmed2./3600;
T2 = [30.7 31.4 31.6 31.9 32.2 32.3 32.5 32.5 32.7 33.3 33.6 33.8 34 34 33.5 33.4 33 33.3 33.3 33.1];
Hr2 = [44.4 40.2 39.5 39.1 39.2 38.6 38.2 38 37.7 36.8 36 35.7 35.7 35.7 36.1 37.1 37.1 36.6 36.8 37];
Vp2 = [3.6 3.61 3.61 3.6 3.59 3.58 3.57 3.56 3.54 3.56 3.56 3.57 3.59 3.59 3.59 3.6 3.6 3.6 3.6 3.6 ];

%Datos totales
%{
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
%}

Ro1 = zeros(1,length(T1));
Roa1 = zeros(1,length(T1));
Mu1 = zeros(1,length(T1));
I1 = zeros(1,length(T1));
pressure1 = zeros(1,length(T1));
eps1 = zeros(1,length(T1));



%% Medicion 1

for i = 1:length(T1)
    %Presión
    
    R = 469; %resistencia [ohm]
    pmax = 2.5; %pmax del rango sensor de presion [bar]
    imin = 4e-3; %corriente minima de salida del sensor [A]
    imax = 20e-3; %corriente maxima de salida del sensor [A]
    
    I1(i) = Vp1(i)/R; %corriente equivalente de presion
    
    pressure1(i) = (I1(i)/imax)*pmax*1e5; %[Pa]
       
    %Densidad
    
    Pamb1 = pressure1(i);
    P1 = 0;
    t1 = T1(i);
    Hum1 = Hr1(i);
    Tamb1 = 26;%[°C]
    
    [ro1,roa1] = density(t1,Hum1,P1,Pamb1); %ro:densidad aire humedo -  roa:densidad aire seco
    [rho1,mu1] = moistAir2(Tamb1,Hum1,P1,Pamb1,t1);%rho:no se usa - mu:viscosidad
 
    Ro1(i) = ro1; %densidad aire húmedo [kg/m3]
    Roa1(i) = roa1; %densidad aire seco [kg/m3]
    Mu1(i) = mu1; %viscosidad
    
    P2_1 = pressure1(i) - dP1(i);
    k = 1.4;
    r1 = P2_1./pressure1(i);
    eps1(i) =  sqrt((r1^(2/k))*(k/(k-1))*((1-r1^((k-1)/k))/(1-r1))*((1-beta^4)/(1-(beta^4)*(r1^(2/k)))));% [dimlss] Factor de correccion por expansion. (Fluidos compresibles)

end

clear i

n = 300;
Cd1 = linspace(0.98,0.995,n);
Qteo1 = zeros(n,length(Qmed1)); %Una fila de Qteo para cada Cd

for i = 1:n
    
    Qteo1(i,:) = Cd1(i) * eps1 * (pi/4) * d^2 .* sqrt(2.*dP1./(Roa1.*(1-beta^4)));%[m3/s] %Caudal volumétrico calculado con los dP medidos
        
end

% plot(dP,Qmed,'bo')

clear i

SL1_1 = zeros(1,n);
sumSL1_1 = zeros(1,length(Qmed1));
SL2_1 = zeros(1,n);
sumSL2_1 = zeros(1,length(Qmed1));

for j = 1:n
    
    for i = 1:length(Qmed1)
        
        sumSL1_1(i) = abs(Qmed1(i)-Qteo1(j,i));
        sumSL2_1(i) = (Qmed1(i)-Qteo1(j,i))^2;

    end
    
    SL1_1(j) = sum(sumSL1_1);
    SL2_1(j) = sum(sumSL2_1);
end

%{
figure
yyaxis left
plot(Cd1,SL1_1,'ro')
ylabel('L1-norm')

hold on

yyaxis right
plot(Cd1,SL2_1,'go')
xlabel('Coeficiente de descarga')
ylabel('L2-norm')

legend('SL1','SL2','Location','North')
title('L-norm vs Cd')
grid minor
%}

Cdfinal1_1 = Cd1(find(SL1_1 == min(SL1_1)))
Cdfinal2_1 = Cd1(find(SL2_1 == min(SL2_1)))

Qteo1_1 = Cdfinal1_1 * eps1 * (pi/4) * d^2 .* sqrt(2.*dP1./(Roa1.*(1-beta^4)));%[m3/s];
Qteo2_1 = Cdfinal2_1 * eps1 * (pi/4) * d^2 .* sqrt(2.*dP1./(Roa1.*(1-beta^4)));%[m3/s];

err1_1 = abs(Qteo1_1 - Qmed1);
err2_1 = abs(Qteo2_1 - Qmed1);


%Caudal másico
Qmteo1_1 = Cdfinal1_1 * eps1 * (pi/4) * d^2 .* sqrt(2.*Roa1.*dP1./((1-beta^4)));%[kg/s];;
Qmteo2_1 = Cdfinal2_1 * eps1 * (pi/4) * d^2 .* sqrt(2.*Roa1.*dP1./((1-beta^4)));%[kg/s];;

QmMed1 = Roa1.*Qmed1;

%plot Qvolumetrico
figure

plot(dP1,Qmed1,'bo-')%,dP,Qteo1,'ro-',dP,Qteo2,'go-')
xlabel('Presión diferencial [Pa]')
ylabel('Caudal volumétrico [m^3/s]')
hold on
% errorbar(dP,Qteo1,err1,'ro-')
errorbar(dP1,Qteo2_1,err2_1,'go-')
% legend('Qmedido','Qteorico1 (L1-norm)','Qteorico2 (L2-norm)','Location','Best')
legend('Q_{medido}','Q_{teórico}','Location','Best')
title({'Caudal volumétrico vs presión diferencial';'\it Medición ascendente'})
grid minor

%plot Qmasico
figure

plot(dP1,QmMed1,'bo-',dP1,Qmteo2_1,'go-')%,dP,Qmteo1,'ro-')
xlabel('Presión diferencial [Pa]')
ylabel('Caudal másico [kg/s]')
hold on
% % errorbar(dP,Qteo1,err1,'ro-')
% errorbar(dP,Qteo2,err2,'go-')
% legend('Qmedido','Qteorico1 (L1-norm)','Qteorico2 (L2-norm)','Location','Best')
legend('Qm_{medido}','Qm_{teórico}','Location','Best')
title({'Caudal másico vs presión diferencial';'\it Medición ascendente'})
grid minor



%
%
n1=2;
ep1_1 = (err1_1./Qmed1)*100;
ep2_1 = (err2_1./Qmed1)*100;
mEp1_1 = mean(ep1_1(n1:end))*ones(1,length(dP1));
mEp2_1 = mean(ep2_1(n1:end))*ones(1,length(dP1));

figure
plot(dP1(n1:end),ep1_1(n1:end),'ro-',dP1(n1:end),ep2_1(n1:end),'go-')
hold on
plot(dP1(n1:end),mEp1_1(n1:end),'r--',dP1(n1:end),mEp2_1(n1:end),'g--')%,'LineWidth',2)
ylabel('\DeltaQ_v [%]')
xlabel('Presión diferencial [Pa]')
title({'\bf Diferencia porcentual entre Q_{teórico} y Q_{medido} vs Presión diferencial';'\it Medición ascendente'})
legend(['L1-norm, C_d = ',num2str(Cdfinal1_1)],['L2-norm, C_d = ',num2str(Cdfinal2_1)],...
    ['\DeltaQ_v promedio = ',num2str(mEp1_1(1)),'%'],['\DeltaQ_v promedio = ',num2str(mEp2_1(1)),'%'])
grid minor

clear i
clear j

%% Medicion 2

Ro2 = zeros(1,length(T2));
Roa2 = zeros(1,length(T2));
Mu2 = zeros(1,length(T2));
I2 = zeros(1,length(T2));
pressure2 = zeros(1,length(T2));
eps2 = zeros(1,length(T2));


for i = 1:length(T2)
    %Presión
    
    R = 469; %resistencia [ohm]
    pmax = 2.5; %pmax del rango sensor de presion [bar]
    imin = 4e-3; %corriente minima de salida del sensor [A]
    imax = 20e-3; %corriente maxima de salida del sensor [A]
    
    I2(i) = Vp2(i)/R; %corriente equivalente de presion
    
    pressure2(i) = (I2(i)/imax)*pmax*1e5; %[Pa]
       
    %Densidad
    
    Pamb2 = pressure2(i);
    P2 = 0;
    t2 = T2(i);
    Hum2 = Hr2(i);
    Tamb2 = 26;%[°C]
    
    [ro2,roa2] = density(t2,Hum2,P2,Pamb2); %ro:densidad aire humedo -  roa:densidad aire seco
    [rho2,mu2] = moistAir2(Tamb2,Hum2,P2,Pamb2,t2);%rho:no se usa - mu:viscosidad
 
    Ro2(i) = ro2; %densidad aire húmedo [kg/m3]
    Roa2(i) = roa2; %densidad aire seco [kg/m3]
    Mu2(i) = mu2; %viscosidad
    
    P2_2 = pressure2(i) - dP2(i);
    k = 1.4;
    r2 = P2_2./pressure2(i);
    eps2(i) =  sqrt((r2^(2/k))*(k/(k-1))*((1-r2^((k-1)/k))/(1-r2))*((1-beta^4)/(1-(beta^4)*(r2^(2/k)))));% [dimlss] Factor de correccion por expansion. (Fluidos compresibles)

end

clear i

n = 300;
Cd2 = linspace(0.97,0.99,n);
Qteo2 = zeros(n,length(Qmed2)); %Una fila de Qteo para cada Cd

for i = 1:n
    
    Qteo2(i,:) = Cd2(i) * eps2 * (pi/4) * d^2 .* sqrt(2.*dP2./(Roa2.*(1-beta^4)));%[m3/s] %Caudal volumétrico calculado con los dP medidos
        
end

% plot(dP,Qmed,'bo')

clear i

SL1_2 = zeros(1,n);
sumSL1_2 = zeros(1,length(Qmed2));
SL2_2 = zeros(1,n);
sumSL2_2 = zeros(1,length(Qmed2));

for j = 1:n
    
    for i = 1:length(Qmed2)
        
        sumSL1_2(i) = abs(Qmed2(i)-Qteo2(j,i));
        sumSL2_2(i) = (Qmed2(i)-Qteo2(j,i))^2;

    end
    
    SL1_2(j) = sum(sumSL1_2);
    SL2_2(j) = sum(sumSL2_2);
end

%{
figure
yyaxis left
plot(Cd2,SL1_2,'ro')
ylabel('L1-norm')

hold on

yyaxis right
plot(Cd2,SL2_2,'go')
xlabel('Coeficiente de descarga')
ylabel('L2-norm')

legend('SL1','SL2','Location','North')
title('L-norm vs Cd')
grid minor
%}

Cdfinal1_2 = Cd2(find(SL1_2 == min(SL1_2)))
Cdfinal2_2 = Cd2(find(SL2_2 == min(SL2_2)))

Qteo1_2 = Cdfinal1_2 * eps2 * (pi/4) * d^2 .* sqrt(2.*dP2./(Roa2.*(1-beta^4)));%[m3/s];
Qteo2_2 = Cdfinal2_2 * eps2 * (pi/4) * d^2 .* sqrt(2.*dP2./(Roa2.*(1-beta^4)));%[m3/s];

err1_2 = abs(Qteo1_2 - Qmed2);
err2_2 = abs(Qteo2_2 - Qmed2);


%Caudal másico
Qmteo1_2 = Cdfinal1_2 * eps2 * (pi/4) * d^2 .* sqrt(2.*Roa2.*dP2./((1-beta^4)));%[kg/s];;
Qmteo2_2 = Cdfinal2_2 * eps2 * (pi/4) * d^2 .* sqrt(2.*Roa2.*dP2./((1-beta^4)));%[kg/s];;

QmMed2 = Roa2.*Qmed2;

%plot Qvolumetrico
figure

plot(dP2,Qmed2,'bo-')%,dP,Qteo1,'ro-',dP,Qteo2,'go-')
xlabel('Presión diferencial [Pa]')
ylabel('Caudal volumétrico [m^3/s]')
hold on
% errorbar(dP,Qteo1,err1,'ro-')
errorbar(dP2,Qteo2_2,err2_2,'go-')
% legend('Qmedido','Qteorico1 (L1-norm)','Qteorico2 (L2-norm)','Location','Best')
legend('Q_{medido}','Q_{teórico}','Location','Best')
title({'Caudal volumétrico vs presión diferencial';'\it Medición descendente'})
grid minor

%plot Qmasico
figure

plot(dP2,QmMed2,'bo-',dP2,Qmteo2_2,'go-')%,dP,Qmteo1,'ro-')
xlabel('Presión diferencial [Pa]')
ylabel('Caudal másico [kg/s]')
hold on
% % errorbar(dP,Qteo1,err1,'ro-')
% errorbar(dP,Qteo2,err2,'go-')
% legend('Qmedido','Qteorico1 (L1-norm)','Qteorico2 (L2-norm)','Location','Best')
legend('Qm_{medido}','Qm_{teórico}','Location','Best')
title({'Caudal másico vs presión diferencial';'\it Medición descendente'})
grid minor


%
%
n2=1;
ep1_2 = (err1_2./Qmed2)*100;
ep2_2 = (err2_2./Qmed2)*100;
mEp1_2 = mean(ep1_2(1:end-n2))*ones(1,length(dP2));
mEp2_2 = mean(ep2_2(1:end-n2))*ones(1,length(dP2));

figure(8)
plot(dP2(1:end-n2),ep1_2(1:end-n2),'ro-',dP2(1:end-n2),ep2_2(1:end-n2),'go-')
hold on
plot(dP2(1:end-n2),mEp1_2(1:end-n2),'r--',dP2(1:end-n2),mEp2_2(1:end-n2),'g--')%,'LineWidth',2)
ylabel('\DeltaQ_v [%]')
xlabel('Presión diferencial [Pa]')
title({'\bf Diferencia porcentual entre Q_{teórico} y Q_{medido} vs Presión diferencial';...
    '\it Medición descendente'})
legend(['L1-norm, C_d = ',num2str(Cdfinal1_2)],['L2-norm, C_d = ',num2str(Cdfinal2_2)],...
    ['\DeltaQ_v promedio = ',num2str(mEp1_2(1)),'%'],['\DeltaQ_v promedio = ',num2str(mEp2_2(1)),'%'])
grid minor

figure
plot(dP2,QmMed2,'bo-',dP2,Qmteo2_2,'go--',dP1,Qmteo1_1,'r>--',dP1,QmMed1,'k>-')
ylabel('Caudal másico [kg/s]')
xlabel('Presión diferencial [Pa]')
legend('Qm medicion descendente','Qm teórico medicion descendente',...
    'Qm teórico medición ascendente','Qm medición ascendente',...
    'Location','SouthEast')
title({'Q_{m} calculado y medido vs \DeltaP';'Medición ascendente y descendente'})
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
%}
