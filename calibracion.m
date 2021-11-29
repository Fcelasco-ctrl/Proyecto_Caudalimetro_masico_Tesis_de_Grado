close all
clear
clc
%{
%% MEDICION FINAL CON SENSORES - 20/12/2018
%Condiciones nominales: Hr = 63%; Tamb = 26°C; Patm = 1009 HPa 
%Medicion ascendente
Rango = [1 1 1 1 2 2 2 2 3 3 3 3 3]; %rango de flujo

Dp = [200 215 240 400 850 1105 1340 1580 1900 2130 2240 2350 2460]; %presion diferencial [Pa]

T = [25.8 27.4 27.7 28.3 30 30.5 31.4 31.8 28.8 30.5 31.2 31.7 32]; %temperatura [°C]

Hr = [69.8 64.4 63.1 61.8 58.8 57.5 55.1 54 65.7 59.8 57.2 55.4 54.3]; %humedad relativa [%]

Vp = [3.706 3.7 3.698 3.697 3.687 3.690 3.693 3.694 3.703 3.705 3.707 3.708 3.71]; %tension para presion absoluta [V]. R = 469 ohm.

Qvmed = [21.64 46.8 53 64 51.5 58 63.8 68.3 38.5 40.3 41.3 42.3 43.1]; %caudal volumetrico medido con el flujometro[m^3/h]

Pmed = [0.8 0.87 0.95 1.44 2.95 3.93 4.95 5.93 6.95 7.95 8.45 8.9 9.4];%presion absoluta medida con el flujómetro [cm H2O]

%Medicion descendente

Dp2 = [207 255 390 700 850 980 1105 1220 1340 1460 1570 1910 2015 2130 2240 2340];

T2 = [29.3 29.5 27.6 31.8 31.7 31.4 30.8 28.6 32.1 31.6 31 30.5 30.2 29.8 29.2 27.8];

Hr2 = [56 57.5 61.6 48.1 48.7 49.3 51.3 54.3 48.4 49.7 50.8 51.9 53 54.1 55.8 60.1];

Vp2 = [3.702 3.710 3.72 3.58 3.58 3.586 3.6 3.592 3.626 3.631 3.635 3.645 3.651 3.655 3.666 3.667];

Qvmed2 = [40.3 46 55.2 40.9 44.5 47.5 50.1 52.5 54.9 57 59 33 34 34.8 35.8 36.4];

Pmed2 = [0.82 0.97 1.45 2.45 2.94 3.45 3.93 4.45 4.92 5.40 5.87 6.94 7.40 7.93 8.40 8.93];



Ro = zeros(1,length(T));
Roa = zeros(1,length(T));
Mu = zeros(1,length(T));
I = zeros(1,length(T));
pressure = zeros(1,length(T));
eps = zeros(1,length(Dp));

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
 
    Ro(i) = ro; %densidad aire húmedo [kg]
    Roa(i) = roa; %densidad aire seco [kg]
    Mu(i) = mu; %viscosidad
    
    
    P2 = pressure(i) - Dp(i);
    k = 1.4;
    beta = 0.4;
    r = P2./pressure(i);
    eps(i) =  sqrt((r^(2/k))*(k/(k-1))*((1-r^((k-1)/k))/(1-r))*((1-beta^4)/(1-(beta^4)*(r^(2/k)))));% [dimlss] Factor de correccion por expansion. (Fluidos compresibles)


   
end

Ro2 = zeros(1,length(T2));
Roa2 = zeros(1,length(T2));
Mu2 = zeros(1,length(T2));
I2 = zeros(1,length(T2));
pressure2 = zeros(1,length(T2));
eps2 = zeros(1,length(Dp2));

for j = 1:length(T2)
    %Presión
    
    R = 469; %resistencia [ohm]
    pmax = 2.5; %pmax del rango sensor de presion [bar]
    imin = 4e-3; %corriente minima de salida del sensor [A]
    imax = 20e-3; %corriente maxima de salida del sensor [A]
    
    I2(j) = Vp2(j)/R; %corriente equivalente de presion
    
    pressure2(j) = (I2(j)/imax)*pmax*1e5; %[Pa]
    
    
    
    %Densidad
    
    Pamb2 = pressure2(j);
    P2 = 0;
    t2 = T2(j);
    Hum2 = Hr2(j);
    Tamb2 = 26;%[°C]
    
    [ro2,roa2] = density(t2,Hum2,P2,Pamb2); %ro:densidad aire humedo -  roa:densidad aire seco
    [rho2,mu2] = moistAir2(Tamb2,Hum2,P2,Pamb2,t2);%rho:no se usa - mu:viscosidad
 
    Ro2(j) = ro2; %densidad aire húmedo [kg/m^3]
    Roa2(j) = roa2; %densidad aire seco [kg/m^3]
    Mu2(j) = mu2; %viscosidad
    
    P22 = pressure2(j) - Dp2(j);
    k = 1.4;
    beta = 0.4;
    r = P22./pressure2(j);
    eps2(j) =  sqrt((r^(2/k))*(k/(k-1))*((1-r^((k-1)/k))/(1-r))*((1-beta^4)/(1-(beta^4)*(r^(2/k)))));% [dimlss] Factor de correccion por expansion. (Fluidos compresibles)

   
end
%}




%Caudal másico
% 
% beta = 0.4;
% D = 3*25.4*1e-3; %[m]
% d = D*beta;
% C = 0.9889; %coeficiente de descarga
% % eps = 0.9926;
% Qm = 0.25*pi*d^2*C.*eps.*sqrt((2.*Ro.*Dp)/(1 - beta^4)); %caudal másico aire humedo [kg/s]
% Qma = 0.25*pi*d^2*C.*eps.*sqrt((2.*Roa.*Dp)/(1 - beta^4)); %caudal másico aire seco [kg/s]
% 
% Qvmedido = Qvmed2./3600; %[m^3/s]  
% Qmmedido = Qvmedido.*Ro2; %[kg/s]


% eps2 = ;
% Qm2 = 0.25*pi*d^2*C.*eps2.*sqrt((2.*Ro2.*Dp2)/(1 - beta^4)); %caudal másico aire humedo [kg/s]
% Qma2 = 0.25*pi*d^2*C.*eps2.*sqrt((2.*Roa2.*Dp2)/(1 - beta^4)); %caudal másico aire seco [kg/s]





% figure
% 
% plot(Qm,Dp,'bo-',Qma,Dp,'ro-')
% xlabel('Caudal másico [kg/s]')
% ylabel('Caída de presión [Pa]')
% 
% title('Caida de presión vs Caudal másico')
% grid minor


%
%% Teórico

Qvc =[ 0.0100    0.0107    0.0114    0.0122    0.0129    0.0136    0.0143    0.0150    0.0157    0.0165    0.0172...
       0.0179    0.0186    0.0193    0.0200    0.0208    0.0215    0.0222    0.0229    0.0236    0.0243    0.0251...
       0.0258    0.0265    0.0272    0.0279    0.0286    0.0294    0.0301    0.0308    0.0315    0.0322    0.0330...
       0.0337    0.0344    0.0351    0.0358    0.0365    0.0373    0.0380    0.0387    0.0394    0.0401    0.0408...
       0.0416    0.0423    0.0430    0.0437    0.0444    0.0451    0.0459    0.0466    0.0473    0.0480    0.0487...
       0.0495    0.0502    0.0509    0.0516    0.0523    0.0530    0.0538    0.0545    0.0552    0.0559    0.0566...
       0.0573    0.0581    0.0588    0.0595    0.0602    0.0609    0.0616    0.0624    0.0631    0.0638    0.0645...
       0.0652    0.0659    0.0667]; %caudal volumetrico teorico

    dPc = 1e3 *[0.1081    0.1242    0.1414    0.1597    0.1791    0.1996    0.2212    0.2440    0.2678    0.2928    0.3189...
                0.3461    0.3744    0.4038    0.4343    0.4660    0.4987    0.5326    0.5676    0.6037    0.6409    0.6792...
                0.7187    0.7592    0.8009    0.8436    0.8875    0.9325    0.9786    1.0259    1.0742    1.1237    1.1742...
                1.2259    1.2787    1.3326    1.3876    1.4437    1.5010    1.5593    1.6188    1.6793    1.7410    1.8038...
                1.8677    1.9328    1.9989    2.0662    2.1345    2.2040    2.2746    2.3463    2.4191    2.4930    2.5681...
                2.6442    2.7215    2.7999    2.8794    2.9600    3.0417    3.1245    3.2085    3.2935    3.3797    3.4670...
                3.5554    3.6449    3.7355    3.8272    3.9201    4.0140    4.1091    4.2053    4.3026    4.4010    4.5005...
                4.6011    4.7029    4.8057]; %delta p teórico


% deltaPCFD = 1.02*[750 1000 2250 2750];
% qCFD = [0.028 0.032 0.048 0.054];
% 
% hold on
% dens = 1.1235;%densidad promedio aire húmedo
% densa = 1.1055;%densidad promedio aire seco
% Qmteor = Qvc.*dens;
% Qmcfd = qCFD*dens;
% Qmateor = Qvc.*densa; %caudal másico teórico aire seco

% plot(Qmteor,dPc,'gs-',Qmcfd,deltaPCFD,'mv-')

% plot(Qmteor,dPc,'gs-',Qmmedido,Dp2,'bo-')
% plot(Qm2,Dp2,'ks-')

% plot(Qmateor,dPc,'rs-')
% legend('Aire húmedo','Aire seco','Qm húmedo teórico','Qm CFD','Aire húmedo 2','Location','Best')


% figure
% plot(Qmteor,dPc,'gs-',Qmmedido,Dp,'bo-')
%}
%
%% MEDICION FINAL CON SENSORES - 21/12/2018
%Condiciones nominales: Hr = 54%; Tamb = 23.8°C; Patm = 1008 HPa 
%Medicion ascendente
Rango = [1 1 1 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3]; %rango de flujo

Dp = [190 240 390 610 714 850 980 1100 1220 1350 1450 1570 1800 1900 2000 2120 2220 2340 2465]; %presion diferencial [Pa]
Dp2=[2450 2350 2220 2130 2020 1900 1890 1680 1570 1460 1338 1225 1100 980 850 710 615 385 248 224];
T = [25.5 26.5 27 27.6 28.3 28.5 28.9 29.3 30 30.4 30.8 30.9 31.2 31.5 31.5 31.7 31.8 31.9 32.1]; %temperatura [°C]
T2 = [30.7 31.4 31.6 31.9 32.2 32.3 32.5 32.5 32.7 33.3 33.6 33.8 34 34 33.5 33.4 33 33.3 33.3 33.1];
Hr = [57 52.2 51.2 50.2 48.9 48.1 47.3 45.9 44.1 42.9 42.3 42.3 41.1 41.4 40.7 39.9 38.4 39.3 39.2]; %humedad relativa [%]
Hr2 = [44.4 40.2 39.5 39.1 39.2 38.6 38.2 38 37.7 36.8 36 35.7 35.7 35.7 36.1 37.1 37.1 36.6 36.8 37];

Vp = [3.61 3.618 3.617 3.615 3.617 3.615 3.61 3.59 3.59 3.6 3.62 3.61 3.60 3.6 3.59 3.58 3.59 3.6 3.59]; %tension para presion absoluta [V]. R = 469 ohm.
Vp = 4*ones(1,length(Dp));
Vp2 = [3.6 3.61 3.61 3.6 3.59 3.58 3.57 3.56 3.54 3.56 3.56 3.57 3.59 3.59 3.59 3.6 3.6 3.6 3.6 3.6 ];
Qvmed = [36.8 53 64.3 81.5 91 99 105.5 111.5 117 122.3 125.5 130.6 144.5 148.3 152.5 156.5 160.2 164 168]; %caudal volumetrico medido con el flujometro[m^3/h]
Qvmed2= [167.5 164 160 156 152 148.5 144.6 140.3 130.3 126.3 121.3 116.3 111.3 105.3 98.5 90.5 81.2 63.5 53 36.8];

Pmed = [0.8 0.98 1.48 2.11 2.43 2.95 3.44 3.93 4.45 4.93 5.4 5.9 6.44 6.95 7.43 7.94 8.92 8.9 9.45];%presion absoluta medida con el flujómetro [cm H2O]
Pmed2 = [9.42 8.92 8.4 7.93 7.41 6.92 6.45 5.94 5.93 5.42 4.92 4.42 3.91 3.45 2.94 2.43 2.16 1.44 0.97 0.9];

Ro = zeros(1,length(T));
Roa = zeros(1,length(T));
Mu = zeros(1,length(T));
I = zeros(1,length(T));
pressure = zeros(1,length(T));
eps = zeros(1,length(Dp));

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
 
    Ro(i) = ro; %densidad aire húmedo [kg]
    Roa(i) = roa; %densidad aire seco [kg]
    Mu(i) = mu; %viscosidad
    
    
    P2 = pressure(i) - Dp(i);
    k = 1.4;
    beta = 0.4;
    r = P2./pressure(i);
    eps(i) =  sqrt((r^(2/k))*(k/(k-1))*((1-r^((k-1)/k))/(1-r))*((1-beta^4)/(1-(beta^4)*(r^(2/k)))));% [dimlss] Factor de correccion por expansion. (Fluidos compresibles)


   
end
% dens = mean(Ro);
% densa = mean(Roa);

beta = 0.4;
D = 3*25.4*1e-3; %[m]
d = D*beta;
C = 0.9889; %coeficiente de descarga

% Qmteor = Qvc.*dens;
Qvmedido = Qvmed./3600; %[m^3/s]
Qvmedido2 = Qvmed2/3600;
Qmmedido = Qvmedido.*Ro; %[kg/s]
Qmmedido2 = Qvmedido2*mean(Ro);
figure
Qmcalculada = 0.25*pi*d^2*C.*eps.*sqrt((2.*Ro.*Dp)/(1 - beta^4)); %caudal másico aire humedo [kg/s]
deltaPCFD = 1.02*[750 1000 2250 2750];
qCFD = [0.028 0.032 0.048 0.054];
yyaxis left
plot(Qmcalculada,Dp,'gs-',Qmmedido,Dp,'bo-',qCFD,deltaPCFD,'mv-'); hold on;
plot(Qmmedido2,Dp2,'b*');
grid minor
title('Caida de presion vs Caudal másico')
xlabel('Caudal másico [kg/s]')
ylabel('Caida de presión [Pa]')
legend('Qm Venturi','Qm Super Flow','Qm CFD','Qm SuperFlow 2','Location','Best')

yyaxis right
diff = (Qmcalculada - Qmmedido);
plot(Qmcalculada,diff,'ro-')

%}
figure
plot(Dp,Qvmedido,'bo')
grid on