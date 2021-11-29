clear all
close all
clc

%% Calculo de factor de fugacidad epsilon:
% Flujeo 1:

dP = [200 215 240 400 850 1105 1340 1580 1900 2130 2240 2350 2460]; %presion diferencial [Pa]
P1 = 1e4*[9.8774 9.8614 9.8561 9.8534 9.8268 9.8348 9.8428 9.8454 9.8694 9.8747 9.8801 9.8827 9.88814]; %presion a la entrada [Pa]
P2 = P1 - dP;
k = 1.4;
beta = 0.4;

r =P2./P1;
e1 = zeros(1,size(dP,2));
for i = 1:size(dP,2)
    
    e1(i) =  sqrt((r(i)^(2/k))*(k/(k-1))*((1-r(i)^((k-1)/k))/(1-r(i)))*((1-beta^4)/(1-(beta^4)*(r(i)^(2/k)))));% [dimlss] Factor de correccion por expansion. (Fluidos compresibles)

end
em1 = mean(e1);

% Flujeo 2: