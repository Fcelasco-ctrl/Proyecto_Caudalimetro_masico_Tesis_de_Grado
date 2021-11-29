close all
% clear
clc
%% Mediciones tomadas con caudalimetro de hilo caliente conectado
%%%%%%%%%%%%%%% Fecha: 07/06/2019
%Presion barometrica = 1020 hPa; Temperatura ambiente = 23.6; Humedad ambiente =
%38%;
% Valores medidos

suma = 33.38 + 46 + 18.5 + 3 + 100.6/2;
%% Medicion Nº1 - Ascendente
%{
Rango = [1 1 2 2 2 2 2 3 3 3 3];
Dp = [227 352 640 865 1085 1295 1515 1760 1972 2175 2390 2515]; %presion diferencial [Pa]

T = [23.6 23.6 23.6 23.6 23.6 23.6 23.6 24.3 24.3 24.3 24.3 24.3]; %temperatura [°C]

Hr = [39 39 39 39 39 39 39 39 39 39 39 39]; %humedad relativa [%]

Vp = [4.94 4.94 4.94 4.94 4.94 4.95 4.95 4.95 4.95 4.96 4.96 4.96]; %tension para presion absoluta [V]. R = 469 ohm.

QMhilo = [0.0141 0.0181 0.0251 0.0290 0.0325 0.0355 0.0382 0.0408 0.0430 0.0451 0.0475 0.0477]; %[kg/s]
%}
%% Medición N°2 - Descendente
%
Dp = [660 985 1400 1865 2295];

QMhilo = [0.0276 0.0337 0.0401 0.0463 0.0518];

Hr = [37 37 37 37 37];

T = [24.5 24.5 24.5 24.5 24.5];

% Vp2 = [];

% P = [102000];
%}
Ro = zeros(1,length(T));
Roa = zeros(1,length(T));
Mu = zeros(1,length(T));
I = zeros(1,length(T));
pressure = zeros(1,length(T));
eps = zeros(1,length(Dp));

pressure = 102000*ones(1,length(Dp));

for i = 1:length(T)
    %Presión
    
    R = 469; %resistencia [ohm]
    pmax = 2.5; %pmax del rango sensor de presion [bar]
    imin = 4e-3; %corriente minima de salida del sensor [A]
    imax = 20e-3; %corriente maxima de salida del sensor [A]
    
    I(i) = Dp(i)/R; %corriente equivalente de presion
    
%     pressure(i) = (I(i)/imax)*pmax*1e5; %[Pa]
    
    
    
    %Densidad
    
    Pamb = pressure(i);
    P = 0;
    t = T(i);
    Hum = Hr(i);
    Tamb = 23.6;%[°C]
    
    [ro,roa] = density(t,Hum,Pamb); %ro:densidad aire humedo -  roa:densidad aire seco
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
C = 0.943; %coeficiente de descarga



Qmcalculada = 0.25*pi*d^2*C.*eps.*sqrt((2.*Ro.*Dp)/(1 - beta^4)); %caudal másico aire humedo [kg/s]
diff = ((Qmcalculada-QMhilo)./QMhilo)*100;

figure
yyaxis left
plot(Dp,Qmcalculada,'gs-',Dp,QMhilo,'bo-')
hold on
% plot(Dp2,QMhilo2,'ks')%,qCFD,deltaPCFD,'mv-'); hold on;

grid on
title('Caudal másico vs Presión diferencial')
ylabel('Caudal másico [kg/s]')
xlabel('Presión diferencial [Pa]')
diff = abs(diff);
yyaxis right
plot(Dp,diff,'r-')
% ylabel('Diferencia porcentual [%]')
legend('Q_m Caudalímetro Venturi','Q_m Caudalímetro Comercial','Location','NorthWest')%,'Qm CFD','Qm SuperFlow 2','Location','Best')
text(Dp(1),diff(1),['\Delta% = ',num2str(round(diff(1),2)),'%']);text(Dp(2)-250,0.85,['\Delta% = ',num2str(round(diff(2),2)),'%']);
text(Dp(3)-250,1.2,['\Delta% = ',num2str(round(diff(3),2)),'%']);text(Dp(4)-250,1.6,['\Delta% = ',num2str(round(diff(4),2)),'%']);
text(Dp(5)-250,1.9,['\Delta% = ',num2str(round(diff(5),2)),'%']);
% annotation('arrow',[0.2 0.15],[0.3 0.1])
% 
% %}
% figure
% plot(Dp,Qvmedido,'bo')
% grid on