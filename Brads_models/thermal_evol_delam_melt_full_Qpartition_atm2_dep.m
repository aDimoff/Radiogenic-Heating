function dTdt = thermal_evol_delam_melt_full_Qpartition_atm2_dep(t,Y,L_m,v_scale_flag,mun,Ts,Rp,rho,g,d,Rc,dT_dP,dcr_crit,DT_sol_dep)

% Solves 14 coupled ordinary differential equations for evolution of mantle
% temperature, crust thickness, heat budget, and CO2 in the atmosphere and
% mantle
% Y(1) is mantle temperature
% Y(2) is the cummulative volume of the crust
% Y(3) is the heat production rate in the crust from U238
% Y(4) is the mantle CO2 reservoir size (in moles)
% Y(5) is the volume of crust present at any one time
% Y(6) is the crust CO2 reservoir size (moles)
% Y(7) is the heat production rate in the mantle from U238
% Y(8) is the thickness of the lithosphere
% Y(9) is the heat production rate in the crust from U235
% Y(10) is the heat production rate in the mantle from U235
% Y(11) is the heat production rate in the crust from Th
% Y(12) is the heat production rate in the mantle from Th
% Y(13) is the heat production rate in the crust from K
% Y(14) is the heat production rate in the mantle from K

dTdt=zeros(14,1);

Ev=300000.0;
Rg=8.314;
alpha=3e-5;
kappa=1e-6;
rho_m=3300.0;
dphi_dP=0.15e-9;
Cp=1250; % J kg^-1 K^-1
k=5; % W m^-1 K^-1
C1=0.5; % heat flow scaling law constant
f=0; % fraction of devolatilized carbon that is trapped and still recycles into mantle

tau_rad2=(-4.46e9/log(1/2))*3600*24*365; % For U238, in seconds
tau_rad235=(-7.04e8/log(1/2))*3600*24*365; % For U235, in seconds
tau_radTh=(-1.4e10/log(1/2))*3600*24*365; % For Th, in seconds
tau_radK=(-1.26e9/log(1/2))*3600*24*365; % For Th, in seconds

% Volume of mantle 
V_man=(4/3)*pi*(Rp^3-Rc^3)-Y(5);
if (Y(5)==0)
    delta_c=0;
else
    delta_c=Rp-(Rp^3-(Y(5)*3)/(4*pi))^(1/3);
end
% Surface area of mantle (excluding crust on top of mantle)
As=4*pi*(Rp-delta_c)^2;

% Surface area of actively convecting mantle
As2=4*pi*(Rp-Y(8))^2;

% Volume of actively convecting mantle 
V_man2=(4/3)*pi*((Rp-Y(8))^3-Rc^3);

% Calculate thermal boundary layer thickness
% Define rayleigh number & theta
Rai=(rho*g*alpha*d^3*(Y(1)-Ts))/(kappa*mun*exp(Ev/(Rg*Y(1))));
theta=(Ev*(Y(1)-Ts))/(Rg*Y(1)^2);

% Calculate heat flux from convection to base of lid 
ql=(C1*k*(Y(1)-Ts)/d)*theta^(-4/3)*Rai^(1/3);

% Temp at base of the lid 
Tl=Y(1)-2.5*(Rg*Y(1)^2)/Ev;

% Concentration of heat producing elements (watts per meter cubed) in crust
% and mantle
if (Y(5)==0)
    x_c=0;
else
    x_c=(Y(3)+Y(9)+Y(11)+Y(13))/Y(5);
end
x_m=(Y(7)+Y(10)+Y(12)+Y(14))/V_man;

% Calculate temperature at base of crust, Tc
Tc=(Ts*(Y(8)-delta_c)+Tl*delta_c)/Y(8)+(x_c*delta_c^2*(Y(8)-delta_c))/(2*k*Y(8))+...
    (x_m/k)*((Y(8)*delta_c)/2 + (delta_c^3)/(2*Y(8))-delta_c^2);

P_lid=Y(8)*rho_m*g;
d_melt=((Y(1)-(1423+DT_sol_dep*(max(1,delta_c/dcr_crit))))/(120e-9-dT_dP))/(rho_m*g);
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
    vol_flux=17.8*3.14*Rp^2.0*vi_dim*phi*(min(d_melt,d_dense)-Y(8))*(1/d);
    %     vol_flux2=64*3.14*Rp^2.0*vi_dim*(d_melt-delta)*(1/d); % Upwelling mantle into melt zone flux
elseif (v_scale_flag == 2)
    % This if for making the upwelling area 1/2 the radius of the
    % cylindrical convection cell
    vi_dim=(kappa/d)*0.05*(Rai/theta)^(2/3);
    %     vi_dim=(kappa/d)*0.05*(Rai*min(1/theta,...
    %         C1^(-3/4)*theta^(1/3)*(theta^(1/3)*(theta-1))^(-1/4)*(Y(2)/(d*Rai))^(1/4)*(d/Y(2))))^(2/3);
    vol_flux=(1/2)*17.8*3.14*Rp^2.0*vi_dim*phi*(min(d_melt,d_dense)-Y(8))*(1/d);
    %     vol_flux2=64*3.14*Rp^2.0*vi_dim*(d_melt-delta)*(1/d); % Upwelling mantle into melt zone flux
else
    vi_dim=(kappa/d)*0.4*(Rai/theta)^(1/2);
    %     vi_dim=(kappa/d)*0.4*(Rai*min(1/theta,...
    %         C1^(-3/4)*theta^(1/3)*(theta^(1/3)*(theta-1))^(-1/4)*(Y(2)/(d*Rai))^(1/4)*(d/Y(2))))^(1/2);
    vol_flux=17.8*3.14*Rp^2.0*vi_dim*phi*(min(d_melt,d_dense)-Y(8))*(1/d);
    %     vol_flux2=64*3.14*Rp^2.0*vi_dim*(d_melt-delta)*(1/d); % Upwelling mantle into melt zone flux
end

frac2=(1-(1-phi)^(1/1e-4));
% For using different distribution coefficients for each HPE
fracU=(1-(1-phi)^(1/0.0012)); % Numbers from Beattie 1993
fracTh=(1-(1-phi)^(1/0.00029)); % From Beattie 1993
fracK=(1-(1-phi)^(1/0.0011)); % From Hart & Brooks 1974, using 60% Ol & 40%  px 

rho_r=2800;% Used as melt density, NOT basalt density in this code
plume_vol_flux=0; %0.1/(1e-9*3600*24*365); % A constant "plume" volume flux of 1 km^3/yr
volc=Y(5);

% This is based on Earth's gravity and is written in terms of depth. Need
% to convert to deal with changing g
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
    z_carb=Y(8)+1000;
end
if (imag(z_carb) ~= 0) 
    z_carb=Y(8)+1000;
end

z_melt=min(d_melt,d_dense); % Making the melt depth equal to depth where melting begins 
% rather than some average, and we will track adiabatic cooling of this
% melt all the way to surface
gamma_ad=2e-8; % K/Pa from Foley & Smye 2018 table
% This is adiabatic gradient of melt subtracted by adiabatic gradient of
% solid mantle, because the mantle adiabat has already been subtracted out
% by using the potential temperature; solid mantle adiabat is what I use in
% the melting law, 3.234e-4 K/m, or 1e-8 K/Pa

% Y(1) is mantle internal temperature
% Real surface temperature stays in melt heat loss, because melt is erupted
% above the seafloor
dTdt(1)=(Y(7)+Y(10)+Y(12)+Y(14))/(V_man2*rho*Cp)-(ql*As2)/(V_man2*rho*Cp)-...
    (vol_flux+plume_vol_flux)*rho_r*(L_m+(Y(1)-Ts-P_melt*gamma_ad)*Cp)/(V_man2*rho*Cp);

% Y(8) is thickness of the stagnant lid, delta
if (Tl>Tc+1)
    dTdt(8)=(1/(rho*Cp*(Y(1)-Tl)))*(-ql-Y(8)*x_m+(k*(Tl-Tc))/(Y(8)-delta_c)+...
        x_m*(Y(8)^2-delta_c^2)/(2*(Y(8)-delta_c)));
else
    dTdt(8)=(1/(rho*Cp*(Y(1)-Tl)))*(-ql-delta_c*x_c*(1/2)+(k*(Tl-Ts))/delta_c);
end

% Y(2) is total volume of crust produced
dTdt(2)=(vol_flux+plume_vol_flux); 

% Y(5) is total volume of crust that is present at any given time 
dTdt(5)=(vol_flux+plume_vol_flux)-...
    (vol_flux+plume_vol_flux-4*pi*(Rp-Y(8))^2*min(0,dTdt(8)))*(tanh((delta_c-Y(8))*20)+1);
if (Y(5)>0 && vol_flux>0)
    % Y(3) is total amount of radiogenic heating in crust from U238 (in Watts)
    dTdt(3)=(Y(7)/V_man)*(vol_flux/phi)*fracU-...
        (Y(3)/(Y(5)+eps))*(vol_flux-4*pi*(Rp-Y(8))^2*min(0,dTdt(8)))*(tanh((delta_c-Y(8))*20)+1)-...
        Y(3)/tau_rad2; 
    % Y(9) is for U235 in crust
    dTdt(9)=(Y(10)/V_man)*(vol_flux/phi)*fracU-...
        (Y(9)/(Y(5)+eps))*(vol_flux-4*pi*(Rp-Y(8))^2*min(0,dTdt(8)))*(tanh((delta_c-Y(8))*20)+1)-...
        Y(9)/tau_rad235; 
    % Y(11) is for Th in crust
    dTdt(11)=(Y(12)/V_man)*(vol_flux/phi)*fracTh-...
        (Y(11)/(Y(5)+eps))*(vol_flux-4*pi*(Rp-Y(8))^2*min(0,dTdt(8)))*(tanh((delta_c-Y(8))*20)+1)-...
        Y(11)/tau_radTh; 
    % Y(13) is for K in crust
    dTdt(13)=(Y(14)/V_man)*(vol_flux/phi)*fracK-...
        (Y(13)/(Y(5)+eps))*(vol_flux-4*pi*(Rp-Y(8))^2*min(0,dTdt(8)))*(tanh((delta_c-Y(8))*20)+1)-...
        Y(13)/tau_radK;
    % Y(7) is total amount of radiogenic heating in mantle from U238 (in Watts)
    dTdt(7)=(Y(3)/(Y(5)+eps))*(vol_flux-4*pi*(Rp-Y(8))^2*min(0,dTdt(8)))*(tanh((delta_c-Y(8))*20)+1)-...
        (Y(7)/V_man)*(vol_flux/phi)*fracU-Y(7)/tau_rad2; 
    % Y(10) is for U235 in mantle
    dTdt(10)=(Y(9)/(Y(5)+eps))*(vol_flux-4*pi*(Rp-Y(8))^2*min(0,dTdt(8)))*(tanh((delta_c-Y(8))*20)+1)-...
        (Y(10)/V_man)*(vol_flux/phi)*fracU-Y(10)/tau_rad235; 
    % Y(12) is for Th in mantle
    dTdt(12)=(Y(11)/(Y(5)+eps))*(vol_flux-4*pi*(Rp-Y(8))^2*min(0,dTdt(8)))*(tanh((delta_c-Y(8))*20)+1)-...
        (Y(12)/V_man)*(vol_flux/phi)*fracTh-Y(12)/tau_radTh; 
    % Y(14) is for K in mantle
    dTdt(14)=(Y(13)/(Y(5)+eps))*(vol_flux-4*pi*(Rp-Y(8))^2*min(0,dTdt(8)))*(tanh((delta_c-Y(8))*20)+1)-...
        (Y(14)/V_man)*(vol_flux/phi)*fracK-Y(14)/tau_radK;     
else
    dTdt(3)=-Y(3)/tau_rad2;
    dTdt(9)=-Y(9)/tau_rad235;
    dTdt(11)=-Y(11)/tau_radTh;
    dTdt(13)=-Y(13)/tau_radK;
    dTdt(7)=-Y(7)/tau_rad2;
    dTdt(10)=-Y(10)/tau_rad235;
    dTdt(12)=-Y(12)/tau_radTh;
    dTdt(14)=-Y(14)/tau_radK;
end
if (z_carb<delta_c && Y(5) > 0 && vol_flux>0)
    % Y(4) is amount of carbon in the mantle
    dTdt(4)=-(Y(4)/V_man)*(vol_flux+plume_vol_flux)*(frac2/phi)+...
    f*(Y(6)/Y(5))*(vol_flux+plume_vol_flux-4*pi*(Rp-Y(8))^2*min(0,dTdt(8)))*(tanh((delta_c-Y(8))*20)+1);
    % Y(6) is amount of carbon in crust
    dTdt(6)=(Y(4)/V_man)*(vol_flux+plume_vol_flux)*(frac2/phi)-...
        f*(Y(6)/Y(5))*(vol_flux+plume_vol_flux-4*pi*(Rp-Y(8))^2*min(0,dTdt(8)))*(tanh((delta_c-Y(8))*20)+1);
elseif (Y(5) > 0 && vol_flux>0)
    dTdt(4)=-(Y(4)/V_man)*(vol_flux+plume_vol_flux)*(frac2/phi)+...
    (Y(6)/Y(5))*(vol_flux+plume_vol_flux-4*pi*(Rp-Y(8))^2*min(0,dTdt(8)))*(tanh((delta_c-Y(8))*20)+1);
    dTdt(6)=(Y(4)/V_man)*(vol_flux+plume_vol_flux)*(frac2/phi)-...
        (Y(6)/Y(5))*(vol_flux+plume_vol_flux-4*pi*(Rp-Y(8))^2*min(0,dTdt(8)))*(tanh((delta_c-Y(8))*20)+1);
else
    dTdt(4)=0;
    dTdt(6)=0;
end
