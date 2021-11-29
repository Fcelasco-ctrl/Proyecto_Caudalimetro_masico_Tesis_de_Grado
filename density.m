%% densidad_error_prop - Densidad del aire humedo con error
%
% Esta funcion devuelve la densidad ro del aire humedo a una dada temperatura t [°C],
% humedad relativa Hum [%] y presion relativa P [Pa]. Pamb es la presion ambiental la cual se utiliza
% para obtener la presion absoluta para realizar el calculo.
%
% densidad
% [ro,err] = densidad(t,Hum,P,Pamb)
% ro [kg/m^3]
% incert [kg/m^3]


function [ro,roa] = density(t,Hum,P)

% Pamb [Pa] presion ambiental
% P [Pa] presion relativa dentro del tubo
% h [%] humedad relativa
% t [°C]
% T [K]

p = P; %[Pa] presion absoluta
T = t + 273.15; %[K]
Ma = 0.028963512440; %[kg/mol] masa molar del aire
Mv = 0.018015; %[kg/mol] masa molar del agua
R = 8.314510; %[J/mol*K] constante universal de los gases
ndatos = length(t);
%% Presion de saturacion de vapor psv

A = 1.2378847e-5; %[1/K^2]
B = -1.9121316e-2; %[1/K]
C = 33.93711047;
D = -6.3431645e3; %[K]

psv = exp(A.*T.^2 + B.*T + C + D./T); %[Pa] presion de saturacion de vapor

%% Fraccion molar del vapor de agua

% Factor de fugacidad f
alfa = 1.00062;
beta = 3.14e-8; %[1/Pa]
gama = 5.6e-7; %[1/K^2]
f = alfa + beta.*p + gama.*t.^2; %factor de fugacidad

h = 0.01*Hum; %humedad relativa

xv = h.*f.*psv./p; %fraccion molar de vapor

%% Factor de compresibilidad Z

a0 = 1.58123e-6*ones(1,ndatos); %[K/Pa]
a1 = -2.9331e-8; %[1/Pa]
a2 = 1.1043e-10; %[1/K*Pa]

b0 = 5.707e-6*ones(1,ndatos); %[K/Pa]
b1 = -2.051e-8;%[1/Pa]

c0 = 1.9898e-4; %[K/Pa]
c1 = -2.376e-6; %[1/Pa]

d = 1.83e-11*ones(1,ndatos); %[K^2/Pa^2]

e = -0.765e-8; %[K^2/Pa^2]

Z = 1 - (p./T).*(a0 + a1.*t + a2.*t.^2 + (b0 + b1.*t).*xv + (c0 + c1.*t).*xv.^2) + ...
    (p./T).^2 .* (d + e.*xv.^2);


ro = (( (p.*Ma)./(Z*R.*T) ) .* (1 - xv.*(1 - Mv/Ma))); %densidad del aire humedo
roa = (1-xv).*( (p.*Ma)./(Z*R.*T) );%densidad del aire seco
%% Error
%{
dT_dt = 1*ones(1,ndatos);
dpsv_dT = (2*A.*T + B*ones(1,ndatos) - D./T.^2).*exp(A*T.^2 + B.*T + C*ones(1,ndatos) + D./T);
df_dp = beta*ones(1,ndatos);
df_dt = 2*gama.*t;
dxv_dh = f.*psv./p;
dxv_df = h.*psv./p;
dxv_dp = -h.*f.*psv./(p.^2);
dxv_dpsv = h.*f./p;
dZ_dp = -(1./T).*(a0+a1.*t+a2.*t.^2+(b0+b1.*t).*xv+(c0+c1.*t).*xv.^2)+(2.*p./(T.^2)).*(d+e.*xv.^2);
dZ_dT = (p./T.^2).*(a0+a1.*t+a2.*t.^2+(b0+b1.*t).*xv+(c0+c1.*t).*xv.^2)-(2.*p.^2./(T.^3)).*(d+e.*xv.^2);
dZ_dt = (-p./T).*(a1+2.*a2.*t+b1.*xv+c1.*xv.^2);
dZ_dxv = (-p./T).*(b0+b1.*t+2*c0.*xv+2*c1.*t.*xv)+(2*p.^2.*e.*xv)./(T.^2);
dro_dp = (Ma./(Z.*R.*T)).*(1-xv.*(1-Mv/Ma));
dro_dZ = (-p.*Ma./(Z.^2.*R.*T)).*(1-xv.*(1-Mv/Ma));
dro_dT = (-p.*Ma./(Z.*R.*T.^2)).*(1-xv.*(1-Mv/Ma));
dro_dxv = ((-p.*Ma)./(Z.*R.*T)).*(1-Mv/Ma);
dro_dR = (-p.*Ma./(Z.*R.^2.*T)).*(1-xv.*(1-Mv/Ma));

%coeficientes de sensibilidad de cada fuente de error

Cp = (dro_dp + dro_dZ.*dZ_dp + dro_dZ.*dZ_dxv.*dxv_df.*df_dp + ...
    dro_dZ.*dZ_dxv.*dxv_dp + dro_dxv.*dxv_df.*df_dp + dro_dxv.*dxv_dp);

Ct = dro_dZ.*dZ_dT.*dT_dt + dro_dZ.*dZ_dt + dro_dZ.*dZ_dxv.*dxv_df.*df_dt + ...
    dro_dZ.*dZ_dxv.*dxv_dpsv.*dpsv_dT.*dT_dt + dro_dT.*dT_dt + ...
    dro_dxv.*dxv_df.*df_dt + dro_dxv.*dxv_dpsv.*dpsv_dT.*dT_dt;

Ch = dro_dZ.*dZ_dxv.*dxv_dh + dro_dxv.*dxv_dh;

CR = dro_dR;

Cec = 1*ones(1,ndatos);

c = [Cp' Ct' Ch' CR' Cec']; %vector de coeficientes de sensibilidad
%}
%% Incertidumbres
%{
%Nivel de confianza y factor de cobertura k
conf_cob = [68.27   1
            90      1.645
            95      1.960
            95.45   2
            99      2.576
            99.73   3];
        
% Presion barometrica

% i- Calibracion del barometro
UB = (0.5/100)*250000; %Incertidumbre expandida (%ERROR A F.S)
k = conf_cob(4,2); %Factor de cobertura
up1 = (UB/k)*ones(1,ndatos-1);

% ii- Resolucion del barometro
dB = 250000/1023; %division de escala [Pa]
%se asume distrib. de prob. rectangular
sigma = sqrt(1/12); %desvio standard
up2 = dB*sigma*ones(1,ndatos-1);

% iii- Variacion de la presion atmosférica durante el periodo de interes
%se asume distrib. de prob. triangular
sigma = sqrt(1/24); 

pi = zeros(1,ndatos-1);
pf = zeros(1,ndatos-1);
for i = 1:ndatos-1
    
    pi(i) = p(i); %presion al inicio de la medicion
    pf(i) = p(i+1); %presion al final de la medicion
%     up3(i) = (pf - pi) .* sigma;
end
up3 = (pf - pi) .* sigma;

% Incertidumbre estandar:

Up = sqrt(up1.^2 + up2.^2 + up3.^2); %[Pa]

% Temperatura ambiente

% i- Calibración del termometro
UT = 0.62/100 ; %incertidumbre expandida (error a fondo de escala del sensor: AOSONG)
k = conf_cob(4,2); %factor de cobertura considerando una distribucion normal
ut1 = (UT/k)*ones(1,ndatos-1);%[°C]

% ii- Resolucion del termometro
dt = 0.1; % [°C] division de escala 
%se asume distribucion de probabilidad rectangular
sigma = sqrt(1/12); %desvio std
ut2 = dt*sigma*ones(1,ndatos-1);%[°C]
% iii- Variacion de la temperatura ambiente durante el periodo de interes

ti = zeros(1,ndatos-1);%[°C] temperatura al inicio de la medicion
tf = zeros(1,ndatos-1);%[°C] temperatura al final de la medicion
for i = 1:ndatos-1
    
    ti(i) = t(i); %temperatura al inicio de la medicion
    tf(i) = t(i+1); %temperatura al final de la medicion
end
%
%se asume una distribucion de probabilidad rectangular
sigma = sqrt(1/24);
ut3 = (tf - ti).*sigma; %[°C]

% Incertidumbre estandar

Ut = sqrt(ut1.^2 + ut2.^2 + ut3.^2); %[°C]

% Humedad relativa del aire

% i- Calibración del higrómetro
UH = (5/100); %Incertidumbre expandida (error a fondo de escala)
%se considera una distribucion normal de probabilidad
k = conf_cob(4,2); %Factor de cobertura
uh1 = (UH/k)*ones(1,ndatos-1); %

% ii- Resolucion del higrometro
dh = 0.1; %resolucion [%] VALIDA PARA SENSOR DE ARDUINO AOSONG
% se considera disribucion triangular de probabilidad
sigma = sqrt(1/12);
uh2 = dh*sigma*ones(1,ndatos-1);
%incertidumbres estandar

%iii- Variacion de la humedad relativa del aire humedo durante el periodo
%de interes
hi = zeros(1,ndatos-1);%humedad al inicio de la medicion
hf = zeros(1,ndatos-1);%humedad relativa al final de la medicion
for i = 1:ndatos-1
    
    hi(i) = h(i); %humedad al inicio de la medicion
    hf(i) = h(i+1); %humedad al final de la medicion
end

%se considera distribucion de probabilidad rectangular
sigma = sqrt(1/24);%desvio std
uh3 = (hf - hi) .* sigma;

%Incertidumbre estandar

Uh = sqrt(uh1.^2 + uh2.^2 + uh3.^2);

% Constante R de los gases

UR = 84e-7*ones(1,ndatos-1); %[J/mol*K]

%Incertidumbre del ajuste de la ecuacion

conf = conf_cob(4,1)/100;
incertid = 1e-4; %[kg/m^3]
Uec = conf*incertid*ones(1,ndatos-1); %[kg/m^3]

u = [Up' Ut' Uh' UR' Uec']; %vector de incertidumbres std


e = zeros(ndatos-1,size(u,2)-2);
for i = 1:ndatos-1
    
    for j = 1:size(u,2)-2
        
        e(i,j) = e(i,j) + abs((c(i,j)*u(i,j))/ro(i))*100;%error porcentual de cada variable en el calculo total
    end
end

rad = zeros(ndatos-1,size(u,2));
for i = 1:ndatos-1
    
    for j = 1:size(u,2)
        
        rad(i,j) = rad(i,j) + (c(i,j)*u(i,j))^2;
        
    end
end



uro = zeros(ndatos-1,1);
for i = 1:ndatos-1
    
    uro(i) = uro(i) + sqrt(sum(rad(i,:)));
end
incert = uro.*2;
error = (incert./ro(1:end-1))*100;
%}
%{
switch variable
    
    case 't'
        
        v = t(1:end-1);
        label = 'Temperatura [°C]';
        
    case 'h'
        
        v = h(1:end-1);
        label = 'Humedad relativa [%]';
        
    case 'p'
        
        v = p(1:end-1);
        label = 'Presión absoluta [Pa]';
end
%}
%{
figure
plot(v,incert)
xlabel(label)
ylabel('Incertidumbre [kg/m^3]')
title('Incertidumbre de la densidad')
grid on
figure
plot(v,error)
xlabel(label)
ylabel('Error porcentual [%]')
title('Error porcentual de la densidad')
grid on
%}
end
