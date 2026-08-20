function [value,isterminal,direction] = myevents_delam_Qpartition_atm2_dep(t,y,Fdegas_crit,mun,v_scale_flag,Ts,Rp,rho,g,d,Rc,dT_dP,dcr_crit,DT_sol_dep)
%
value(1:4)=0;

Ev=300000.0;
Rg=8.314;
alpha=3e-5;
kappa=1e-6;
rho_m=3300.0;
dphi_dP=0.15e-9;
% Cp=1250; % J kg^-1 K^-1
% Rc=3488.1*1000; % Core radius
k=5; % W m^-1 K^-1
% Supply limited weathering rate equation
% First need moles of CaO, FeO, MgO, K2O, Na2O per 100 g 
% Weight percents: MgO 7.58; CaO 11.39; FeO 10.43; Na2O 2.79; K2O 0.16
% Mol/100 g: MgO 0.19; CaO 0.2; FeO 0.145; Na2O 0.045; K2O 0.0017
weath_demand=10*(0.19+0.2+0.145+0.045+0.0017);
A_planet=4*pi*Rp^2;
A_earth=4*pi*6378100.0^2;
f=0;

% Volume of mantle 
V_man=(4/3)*pi*(Rp^3-Rc^3)-y(5);
if (y(5)==0)
    delta_c=0;
else
    delta_c=Rp-(Rp^3-(y(5)*3)/(4*pi))^(1/3);
end

% Calculate thermal boundary layer thickness
% Define rayleigh number & theta
Rai=(rho*g*alpha*d^3*(y(1)-Ts))/(kappa*mun*exp(Ev/(Rg*y(1))));
theta=(Ev*(y(1)-Ts))/(Rg*y(1)^2);

% Temp at base of the lid 
Tl=y(1)-2.5*(Rg*y(1)^2)/Ev;

% Concentration of heat producing elements (watts per meter cubed) in crust
% and mantle
if (y(5)==0)
    x_c=0;
else
    x_c=(y(3)+y(9)+y(11)+y(13))/y(5);
end
x_m=(y(7)+y(10)+y(12)+y(14))/V_man;

% Calculate temperature at base of crust, Tc
Tc=(Ts*(y(8)-delta_c)+Tl*delta_c)/y(8)+(x_c*delta_c^2*(y(8)-delta_c))/(2*k*y(8))+...
    (x_m/k)*((y(8)*delta_c)/2 + (delta_c^3)/(2*y(8))-delta_c^2);

P_lid=y(8)*rho_m*g;
d_melt=((y(1)-(1423+DT_sol_dep*(max(1,delta_c/dcr_crit))))/(120e-9-dT_dP))/(rho_m*g);
P_melt1=rho_m*g*d_melt;
P_melt=min(P_melt1,10e9); % cutoff on melting depth, melt dense at ~ 10 GPa
d_dense=10e9/(rho_m*g);
if (P_melt > P_lid)
    phi=0.5*(P_melt1-P_lid)*dphi_dP;
else
    phi=0;
end
if (v_scale_flag == 1)
        vi_dim=(kappa/d)*0.05*(Rai/theta)^(2/3);
%     vi_dim=(kappa/d)*0.05*(Rai*min(1/theta,...
%         C1^(-3/4)*theta^(1/3)*(theta^(1/3)*(theta-1))^(-1/4)*(Y(2)/(d*Rai))^(1/4)*(d/Y(2))))^(2/3);
    vol_flux=17.8*3.14*Rp^2.0*vi_dim*phi*(min(d_melt,d_dense)-y(8))*(1/d);
%     vol_flux2=64*3.14*Rp^2.0*vi_dim*(d_melt-delta)*(1/d); % Upwelling mantle into melt zone flux
elseif (v_scale_flag == 2)
    % This if for making the upwelling area 1/2 the radius of the
    % cylindrical convection cell
    vi_dim=(kappa/d)*0.05*(Rai/theta)^(2/3);
    %     vi_dim=(kappa/d)*0.05*(Rai*min(1/theta,...
    %         C1^(-3/4)*theta^(1/3)*(theta^(1/3)*(theta-1))^(-1/4)*(Y(2)/(d*Rai))^(1/4)*(d/Y(2))))^(2/3);
    vol_flux=(1/2)*17.8*3.14*Rp^2.0*vi_dim*phi*(min(d_melt,d_dense)-y(8))*(1/d);
    %     vol_flux2=64*3.14*Rp^2.0*vi_dim*(d_melt-delta)*(1/d); % Upwelling mantle into melt zone flux
else
    vi_dim=(kappa/d)*0.4*(Rai/theta)^(1/2);
%     vi_dim=(kappa/d)*0.4*(Rai*min(1/theta,...
%         C1^(-3/4)*theta^(1/3)*(theta^(1/3)*(theta-1))^(-1/4)*(Y(2)/(d*Rai))^(1/4)*(d/Y(2))))^(1/2);
    vol_flux=17.8*3.14*Rp^2.0*vi_dim*phi*(min(d_melt,d_dense)-y(8))*(1/d);
%     vol_flux2=64*3.14*Rp^2.0*vi_dim*(d_melt-delta)*(1/d); % Upwelling mantle into melt zone flux
end

% For assuming that some melt stops at base of the lid or otherwise doesnt
% contribute to mantle cooling, degassing, or crustal growth
%vol_flux=0.5*vol_flux;

volc=y(5);

A=3.125e-3/(rho_m*9.8);
B=835.5;
%crust decarbonation degassing flux
if (delta_c > 0)
    % Just calculate it using my temperature profile-decarbonation temp
    % solution
    % If z_carb > delta_c, then no decarbonation is occuring, and use
    % mantle recycling formulation
    z_carb=(delta_c/2)+(k*(Tc-Ts))/(delta_c*x_c)-(A*rho_m*g*k)/x_c - ...
        (k/x_c)*sqrt(((x_c*delta_c)/(2*k)+(Tc-Ts)/delta_c-A*rho_m*g)^2 + (2*x_c)/k*(Ts-B));
else
    z_carb=y(8)+1000;
end
if (imag(z_carb) ~= 0) 
    z_carb=y(8)+1000;
end
% z_carb=y(8)+1000;
if (delta_c < z_carb)
    volc=y(5);
else
    volc=(4/3)*pi*(Rp^3-(Rp-z_carb)^3);
end

if (y(5) > 0)
    F_dcarb=(1-f)*(y(6)/volc)*((vol_flux)*(1/2))*(tanh((delta_c-z_carb)*20)+1);
    F_dcarb1=(1-f)*(y(6)/volc)*(1/2)*(tanh((delta_c-z_carb)*20)+1);
else
    F_dcarb=0;
    F_dcarb1=0;
end

% 12/1/16: This one needs to have the factor of 1/2, because we do have
% crust thicker than z_carb (often!)

frac2=(1-(1-phi)^(1/1e-4));

if (phi>0)
    Fdegas=(frac2/phi)*vol_flux*y(4)/V_man;
    Fdegas1=(frac2/phi)*y(4)/V_man;
else
    Fdegas=0;
    Fdegas1=0;
end
% Use rho_m for simplicity; we are just using mantle density (at surface
% pressure) for everything related to crust
Fws1=0.1*rho_m*weath_demand;
Fws=0.1*vol_flux*rho_m*weath_demand;

value(1) = Fdegas+F_dcarb-Fdegas_crit*(A_planet/A_earth); % time when degassing first falls below modern Earth rate
value(2) = Fdegas+F_dcarb-0.1*Fdegas_crit*(A_planet/A_earth); % time when degassing first falls below 10 % modern Earth rate
value(3) = Fdegas+F_dcarb; % time when degassing ends
value(4) = Fdegas+F_dcarb-Fws-1e14/(3600*24*365); % Time when weathering becomes supply limited 
value(5) = Fdegas+F_dcarb-0.1*Fdegas_crit*(A_planet/A_earth); % time when degassing first goes above 10 % modern Earth rate
isterminal(1:2) = 0; % stop the integration
isterminal(4)=0;
isterminal(5)=0;
isterminal(3)=0;
direction(1:3) = -1;
direction(4)=1;
direction(5)=1;


