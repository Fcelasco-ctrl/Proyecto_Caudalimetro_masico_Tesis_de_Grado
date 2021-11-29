%% Temp [°C] Hum [%]  P [Pa] Pamb[Pa] T[°C] ndatos [cantidad de datos]
%% P=0 y T=Tamb


function [ rho,mu ] = moistAir2(Tamb,Hum,P,Pamb,T)


Pvs=10.^(7.5.*(Tamb+273.15-273.16)./(Tamb+273.15-35.85)+2.7858); %Presión de saturación de vapor de agua [Pa]
Pv=0.01.*Hum.*Pvs;           %Presión de vapor de agua en ensayo
x=0.622.*Pv./(P+Pamb-Pv);         %Humedad absoluta
rho=(1+x)./(461.56*(0.62198+x)).*(P+Pamb)./(T+237.15);

T=T+273.15;
mua=(0.40401+0.074582*T-5.7171e-5*T.^2+2.9928e-8*T.^3-6.2524e-12*T.^4);
muv=((T/647.27).^0.5)./(0.0181583+0.0177624*(647.27./T)+0.0105287*(647.27./T).^2-0.0036744*(647.27./T).^3);
phiav=(1+(mua./muv).^0.5*(18/29)^0.25).^2/(2^(3/2)+(1+29/18)^0.5);
phiva=(1+(muv./mua).^0.5*(29/18)^0.25).^2/(2^(3/2)+(1+18/29)^0.5);
mu=(mua./(1+phiav.*x*1.61)+muv./(1+phiva./(x*1.61)))*1e-6;



end

