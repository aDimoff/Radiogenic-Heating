close all
clear variables

M_E=5.97e24; % Mass of the Earth in kg 

mass=1.000; % Planet mass (in Earth masses) (parameter to change)
Rp=6378100.0*mass^0.27; % Radius in m
rho=4450.0*mass^0.2; % Mantle density in kg/m^3
g=mass*M_E*6.67e-11/Rp^2; % Surface gravity in m/s^2
d=2890000*mass^0.28; % Mantle thickness in m
Rc=Rp-d; % Core radius

Ev=300000.0; % Vicsosity activation energy in J/mol
Rg=8.314; % Universal gas constant
alpha=3e-5; % Thermal expansion coefficient
kappa=1e-6; % Thermal diffusivity 
rho_m=3300.0; % Mantle density at 0 pressure
dphi_dP=0.15e-9; % Change in melt fraction for a change in melting pressure
Cp=1250; % J kg^-1 K^-1
k=5; % W m^-1 K^-1
C1=0.5; % heat flow scaling law constant
Ts=273.0; % Surface temperature
L_m=600*1e3; % Latent heat in J/kg
dT_dP=2e-8; % adiabatic gradient in mantle in K/Pa
gamma_ad=2e-8; % K/Pa from Foley & Smye 2018 table

f=0;

Ts_star=285; % Present day surface temperature
T_bot_star=1.02*Ts_star-16.7; % Reference ocean bottom temperature
pCO2_star=30; % Present day atmospheric CO2, in Pa
A_planet=4*pi*Rp^2;
E_sfw=75000; % J/mol, activation energy for seafloor weathering
F_sfw_star=0.5e12; % mol of C/yr; present day Earth seafloor weathering flux
F_sfw_star=F_sfw_star/(3600*24*365); % Converted to s
vol_flux_star=23*1e9; % 23 km^3/yr present day Earth melt flux at ridges
vol_flux_star=vol_flux_star/(3600*24*365); % Converted to s
vex=0.5; % exponent for dependence of volanism rate on seafloor weathering

V_man=(4/3)*pi*(Rp^3-Rc^3); % Volume of the mantle
V_man_Earth=(4/3)*pi*(6378100^3-3488100^3); % Volume of Earth's mantle
m_man=V_man*rho; % Mass of the mantle 
m_man_Earth=3.9e24; % mass of Earth's mantle in kg

% Target CMF
target_CMF=0.325;   % Example: Earth-like

% Solve for new core radius Rc_new
% CMF = m_core / M_total = (M_total - m_man) / M_total
% m_man = rho * (4/3)*pi*(Rp^3 - Rc^3)
% So: target_CMF = 1 - [rho * (4/3)*pi*(Rp^3 - Rc^3)] / M_total

% Rearranging to solve for Rc_new:
Rc_new=(Rp^3-(mass*(1-target_CMF))/(rho*(4/3)*pi))^(1/3);

% Update your variables
Rc=Rc_new;
d=Rp-Rc;                     % New mantle thickness
V_man=(4/3)*pi*(Rp^3-Rc^3);  % New mantle volume
m_man=V_man*rho;             % New mantle mass
fprintf('New mantle mass = %.3e kg\n', m_man);
m_core=mass-m_man;        % New core mass
fprintf('New core mass = %.3e kg\n', m_core);
CMF=m_core/mass;          % Should equal target_CMF

fprintf('Adjusted to CMF = %.3f\n', CMF);
fprintf('New core radius = %.1f km\n', Rc/1000);
fprintf('New mantle thickness = %.1f km\n', d/1000);

% % Determine starting radiogenic heating rates for each radionuclide
Earth_U238=0.022*0.9927; %ppm
Earth_U235=0.022*0.0072; %ppm
Earth_Th=0.083; %ppm
Earth_K=261*(0.0117/100); %ppm

Q_U238=9.17e-5; % W/kg of U238
Q_U235=5.75e-4; % W/kg of U235
Q_Th=2.56e-5; % W/kg of Th
Q_K=2.97e-5; % W/kg of K

% Half lives and decay constants for four major HPEs
tau_half_U238=4.46e9; % years
lambda_U238=-tau_half_U238/log(1/2);
tau_half_U235=7.04e8; % years
lambda_U235=-tau_half_U235/log(1/2);
tau_half_Th=1.4e10; % years
lambda_Th=-tau_half_Th/log(1/2);
tau_half_K=1.26e9; % years
lambda_K=-tau_half_K/log(1/2);

% Parameters for lognormal distributions of each HPE based on observations
% in stars
mu_U238=-0.072334;
sigma_U238=0.38916;
mu_U235=-0.072334;
sigma_U235=0.38916;
mu_Th=0.18834;
sigma_Th=0.25947;
mu_K=0.14086;
sigma_K=0.60956;

t_end=10e9*3600*24*365; % End time of the model

Ae=4*pi*Rp^2; % area of planet
Ae0=4*pi*6378100^2; % area of Earth
m_bar_co2=44/1000; % molar mass of CO2
m_bar_w=18/1000; % molar mass of H20
plume_vol_flux=0;
rho_r=2800;% Used as melt density, NOT basalt density in this code
Fdstar=6e12/(3600*24*365);

dcr_crit=0.2*V_man/Ae;
DT_sol_dep=150; % Increase in solidus temperature due to full depletion

% Supply limited weathering rate equation
% First need moles of CaO, FeO, MgO, K2O, Na2O per 100 g 
% Weight percents: MgO 7.58; CaO 11.39; FeO 10.43; Na2O 2.79; K2O 0.16
% Weight percents from Gale et al, 2013
% Mol/100 g: MgO 0.19; CaO 0.2; FeO 0.145; Na2O 0.045; K2O 0.0017
weath_demand=10*(0.19+0.2+0.145+0.045+0.0017); % mol/kg of CO2 drawdown

ex_time=1;
rng('shuffle');

% % Sample from uniform distributions of ref viscosity, C budget, initial T
% mun=10^unifrnd(9.3444,13.3444); % sample from mun=4e10 to 4e12;
% C_tot=10^unifrnd(20,23.5)*(m_man/m_man_Earth);
% Ti_init=unifrnd(1700,2000);
% mu_ref=mun*exp(Ev/(Rg*1623));

% Or can set initial mantle temperature, reference viscosity, and mantle C
% budget directly
% Parameters to vary
Ti_init=2000; % Initial mantle temperature (K)
mu_ref=1e21; % Mantle reference viscosity (Pa s)
mun=mu_ref/exp(Ev/(Rg*1623));
% Total carbon budget of the mantle, scaled by mantle mass relative to
% Earth (in moles)
C_tot=1e22*(m_man/m_man_Earth);

% % Randomly sample from log normal distributions of HPEs
% U238=lognrnd(mu_U238,sigma_U238);
% U235=lognrnd(mu_U235,sigma_U235);
% Th=lognrnd(mu_Th,sigma_Th);
% K=lognrnd(mu_K,sigma_K);

% Or set them to Earth's value (this would be == 1.0)
% from Roederer's paper:
% max enhanced value [Eu/Fe] = +2.45 dex = 281...?
% max enhanced value [U/Fe] = +2.00 = 100
% max enhanced value [Th/Fe] = +2.40 = 251
% enhanced in [K/Fe] = +0.15 = 1.41
% Parameter to vary
U238=1.000;
U235=1.000;
Th=1.000;
K=1.000;

% Determine initial heating rate from each HPE
Q0_U238=m_man*Q_U238*U238*Earth_U238*(1/1e6)*exp(4.5e9/lambda_U238);
Q0_U235=m_man*Q_U235*U235*Earth_U235*(1/1e6)*exp(4.5e9/lambda_U235);
Q0_Th=m_man*Q_Th*Th*Earth_Th*(1/1e6)*exp(4.5e9/lambda_Th);
Q0_K=m_man*Q_K*K*Earth_K*(1/1e6)*exp(4.5e9/lambda_K);
Q0_tot=Q0_U238+Q0_U235+Q0_Th+Q0_K;

Q0_U238_save=Q0_U238;
Q0_U235_save=Q0_U235;
Q0_Th_save=Q0_Th;
Q0_K_save=Q0_K;
Q0_tot_save=Q0_tot;

% Find initial lithosphere thickness with crust thickness of zero
Rai=(rho*g*alpha*d^3*(Ti_init-Ts))/(kappa*mun*exp(Ev/(Rg*Ti_init)));
theta=(Ev*(Ti_init-Ts))/(Rg*Ti_init^2);
V_man=(4/3)*pi*(Rp^3-Rc^3);
ql=(C1*k*(Ti_init-Ts)/d)*theta^(-4/3)*Rai^(1/3);
Tl=Ti_init-2.5*(Rg*Ti_init^2)/Ev;
x_m=Q0_tot/V_man;
func1=@(x) Tl-Ts-ql*x/k - x_m*x^2/k + x_m*x^2/(2*k);
delta_init=fzero(func1,100000);

options = odeset('RelTol',1e-11,'NonNegative',[1 2 3 4 5 6 7 8 9 10 11 12 13 14],...
    'Events',@(t,y) myevents_delam_Qpartition_atm2_dep(t,y,Fdstar,mun,1,Ts,Rp,rho,g,d,Rc,dT_dP,dcr_crit,DT_sol_dep));
% options3 = odeset('RelTol',1e-10,'AbsTol',[1e-10,1e-10,1e5],'NonNegative',3);
options2 = optimset('TolFun',1e-12);

tic
% Use ODE solver to solve model
[T,Y,Te,Ye,ie]=ode15s(@(t,y) thermal_evol_delam_melt_full_Qpartition_atm2_dep(t,y,L_m,1,mun,Ts,Rp,rho,g,...
    d,Rc,dT_dP,dcr_crit,DT_sol_dep),[0 t_end],[Ti_init 0 0 1.0*C_tot 0 0.0*C_tot Q0_U238 ...
    delta_init 0 Q0_U235 0 Q0_Th 0 Q0_K],options);
toc
time(1:length(T))=T(1:length(T))/(1e9*3600*24*365); % Time in Gyrs
Ti(1:length(T))=Y(1:length(T),1); % Mantle temperature 
crust(1:length(T))=Y(1:length(T),2); % Cummulative volume of crust produced
crust2(1:length(T))=Y(1:length(T),5); % Volume of crust present as a function of time
Rman(1:length(T))=Y(1:length(T),4); % Total amount of CO2 (in moles) in mantle
Rman2=Rman;
Rcrust(1:length(T))=Y(1:length(T),6); % Total amount of CO2 (in moles) in the crust
Rcrust2=Rcrust;
delta(1:length(T))=Y(1:length(T),8); % Lithosphere thickness in meters
u238_crust(1:length(T))=Y(1:length(T),3); % Heat production (in Watts) from U238 in the crust
u235_crust(1:length(T))=Y(1:length(T),9); % Heat production (in Watts) from U235 in the crust
th_crust(1:length(T))=Y(1:length(T),11); % Heat production (in Watts) from Th in the crust
K_crust(1:length(T))=Y(1:length(T),13); % Heat production (in Watts) from K in the crust
u238_man(1:length(T))=Y(1:length(T),7); % Heat production (in Watts) from U238 in the mantle
u235_man(1:length(T))=Y(1:length(T),10); % Heat production (in Watts) from U235 in the mantle
th_man(1:length(T))=Y(1:length(T),12); % Heat production (in Watts) from Th in the mantle
K_man(1:length(T))=Y(1:length(T),14); % Heat production (in Watts) from K in the mantle

% Sorting events
supply_lim_start=10;
fdegas1=0;
fdegas2=0;
fdegas_end=0;
for i=1:1:length(Te)
    if ie(i) == 1
        fdegas1=Te(i)/(3600*24*365*1e9);
    elseif ie(i) == 2
        fdegas2=Te(i)/(3600*24*365*1e9);
    elseif ie(i) == 3
        fdegas_end=Te(i)/(3600*24*365*1e9);
    elseif ie(i) == 4
        supply_lim_start=Te(i)/(3600*24*365*1e9);
    end
end

for l=1:1:length(time)
    
    % Volume of mantle
    V_man=(4/3)*pi*(Rp^3-Rc^3)-crust2(l);
    if (crust2(l)==0)
        delta_c=0;
    else
        delta_c=Rp-(Rp^3-(crust2(l)*3)/(4*pi))^(1/3);
    end
    
    % Calculate thermal boundary layer thickness
    % Define rayleigh number & theta
    Rai=(rho*g*alpha*d^3*(Ti(l)-Ts))/(kappa*mun*exp(Ev/(Rg*Ti(l))));
    theta=(Ev*(Ti(l)-Ts))/(Rg*Ti(l)^2);
    P_lid=delta(l)*rho_m*g;
    d_melt=((Ti(l)-(1423+DT_sol_dep*(max(1,delta_c/dcr_crit))))/(120e-9-dT_dP))/(rho_m*g);
    P_melt1=rho_m*g*d_melt;
    P_melt=min(P_melt1,10e9); % cutoff on melting depth, melt dense at ~ 10 GPa
    d_dense=10e9/(rho_m*g);
    d_melt_plot(l)=min(d_melt,d_dense);
    if (P_melt > P_lid)
        phi=0.5*(P_melt1-P_lid)*dphi_dP;
    else
        phi=0;
    end
    vi_dim=(kappa/d)*0.05*(Rai/theta)^(2/3);
    vi_dim_plot(l)=vi_dim;
    %     vi_dim=(kappa/d)*0.05*(Rai*min(1/theta,...
    %         C1^(-3/4)*theta^(1/3)*(theta^(1/3)*(theta-1))^(-1/4)*(Y(2)/(d*Rai))^(1/4)*(d/Y(2))))^(2/3);
    vol_flux=17.8*3.14*Rp^2.0*vi_dim*phi*(min(d_melt,d_dense)-delta(l))*(1/d);
    %     vol_flux2=64*3.14*Rp^2.0*vi_dim*(d_melt-delta)*(1/d); % Upwelling mantle into melt zone flux
    vol_flux_plot(l)=vol_flux;
    
    q_magma(l)=vol_flux*rho_r*(L_m+(Ti(l)-Ts-P_melt*gamma_ad)*Cp)/A_planet;
    ql(l)=(C1*k*(Ti(l)-Ts)/d)*theta^(-4/3)*Rai^(1/3);

    % For partitioning of CO2 into melt
    frac2=(1-(1-phi)^(1/1e-4));
    
    % Temp at base of the lid
    Tl=Ti(l)-2.5*(Rg*Ti(l)^2)/Ev;
    
    Q_crust_tot=u238_crust(l)+u235_crust(l)+th_crust(l)+K_crust(l);
    Q_crust_tot_plot(l)=Q_crust_tot;
    Q_man_tot=u238_man(l)+u235_man(l)+th_man(l)+K_man(l);
    Q_man_tot_plot(l)=Q_man_tot;
    % Concentration of heat producing elements (watts per meter cubed) in crust
    % and mantle
    if (crust2(l)==0)
        x_c=0;
    else
        x_c=Q_crust_tot/crust2(l);
    end
    x_m=Q_man_tot/V_man;
    x_m_plot(l)=x_m;
    x_c_plot(l)=x_c;
    
    % Calculate temperature at base of crust, Tc
    Tc=(Ts*(delta(l)-delta_c)+Tl*delta_c)/delta(l)+(x_c*delta_c^2*(delta(l)-delta_c))/(2*k*delta(l))+...
        (x_m/k)*((delta(l)*delta_c)/2 + (delta_c^3)/(2*delta(l))-delta_c^2);
    
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
        z_carb=delta(l)+1000;
    end
    if (imag(z_carb) ~= 0)
        z_carb=delta(l)+1000;
    end
    if (delta_c < z_carb)
        volc=crust2(l);
    else
        volc=(4/3)*pi*(Rp^3-(Rp-z_carb)^3);
    end

    z_carb_plot(l)=z_carb;
    
    if time(l) > supply_lim_start
%         Rman2(l)=Rman2(l-1)-(Fdegas(l-1))*3600*24*365*1e9*(time(l)-time(l-1));
%         Rcrust2(l)=Rcrust2(l-1)-(F_dcarb(l-1)-Fws1(l-1))*3600*24*365*1e9*(time(l)-time(l-1));
%         pCO2(l)=pCO2(l-1)+(F_dcarb(l-1)+Fdegas(l-1)-Fws1(l-1))*3600*24*365*1e9*(time(l)-time(l-1))*m_bar_co2*g/Ae;
%         Ts_calc(l)=Ts_star+5.6*(log(pCO2(l)/pCO2_star)/log(2));
        if (crust2(l) > 0)
            F_dcarb(l)=(1-f)*(Rcrust2(l)/volc)*((vol_flux)*(1/2))*(tanh((delta_c-z_carb)*20)+1);
        else
            F_dcarb(l)=0;
        end
        if (phi>0)
            Fdegas(l)=vol_flux*(Rman2(l)/V_man)*(frac2/phi);
        else
            Fdegas(l)=0;
        end
        % Use rho_m for simplicity; we are just using mantle density (at surface
        % pressure) for everything related to crust
        Fws1(l)=(0.1*vol_flux*rho_m*weath_demand);
    else
        if (crust2(l) > 0)
            F_dcarb(l)=(1-f)*(Rcrust(l)/volc)*((vol_flux)*(1/2))*(tanh((delta_c-z_carb)*20)+1);
        else
            F_dcarb(l)=0;
        end
        if (phi>0)
            Fdegas(l)=vol_flux*(Rman(l)/V_man)*(frac2/phi);
        else
            Fdegas(l)=0;
        end
        % Use rho_m for simplicity; we are just using mantle density (at surface
        % pressure) for everything related to crust
        Fws1(l)=(0.1*vol_flux*rho_m*weath_demand);
    end
    
end
if max(Fdegas+F_dcarb-Fws1)*3600*24*365>1e14
    supply_lim2=1;
else
    supply_lim2=0;
end

% % Can use this to write results out to a file
% fid=fopen('Earth_example_Ctot1e23','w');
% for l=1:1:length(time)
%     fprintf(fid,'%8.6e %8.6e %8.6e %8.6e %8.6e %8.6e %8.6e\r\n',...
%         time(l),Fdegas(l)*3600*24*365,F_dcarb(l)*3600*24*365,(Fdegas(l)+F_dcarb(l))*3600*24*365,...
%         Fws1(l)*3600*24*365), Ts_calc(l),pCO2(l));
% end
% fclose(fid);

fid=fopen('CMF_sweep/Earth_mantle_temp_CMF_0.325mmmmmmm.dat','w');
for l=1:1:length(time)
    fprintf(fid,'%8.6e %8.6e\r\n',...
        time(l),Ti(l));
end
fclose(fid);

figure(1)
plot(time,Ti)
xlabel('Time [Gyrs]');
ylabel('Mantle temperature [K]');

figure(2)
semilogy(time,Fdegas*3600*24*365,time,Fws1*3600*24*365,time,F_dcarb*3600*24*365,...
    time,(Fdegas+F_dcarb)*3600*24*365,fdegas1,6e12*(Ae/Ae0),'s',fdegas2,0.1*6e12*(Ae/Ae0),...
    'd',fdegas_end,1e6,'^',supply_lim_start,1e15,'o')
legend('Mantle outgassing','Supply limit','Metamorphic outgassing','Total outgassing rate','location','best')
xlabel('Time [Gyrs]');
ylabel('Outgassing rate [mol/Myr]');
% 
% figure(3)
% plot(time,delta/1000)
% xlabel('Time [Gyrs]');
% ylabel('Lithosphere thickness [m]');
% 
% figure(4)
% plot(time,Q_crust_tot_plot/1e12,time,Q_man_tot_plot/1e12)
% xlabel('Time [Gyrs]');
% ylabel('Heat production rate [TW]');
% legend('Crust','Mantle','location','best')

figure(5)
semilogy(time,Rman,time,Rcrust)
xlabel('Time [Gyrs]');
ylabel('Carbon reservoir size [mol]');
legend('Mantle','Crust','location','best')

fid=fopen('CMF_sweep/Earth_lith_crust_melt_CMF_0.325mmmmmmm.dat','w');
for l=1:1:length(time)
    fprintf(fid,'%8.6e %8.6e %8.6e %8.6e %8.6e\r\n',...
        time(l),delta(l)/1e3,(Rp-(Rp^3-(crust2(l)*3)/(4*pi)).^(1/3))/1e3,d_melt_plot(l)/1e3,z_carb_plot(l)/1e3);
end
fclose(fid);

figure(6)
plot(time,delta/1e3,time,(Rp-(Rp^3-(crust2*3)/(4*pi)).^(1/3))/1e3,time,d_melt_plot/1e3,time,z_carb_plot/1e3)
xlabel('Time [Gyrs]');
ylabel('Depth [km]');
legend('Base of lithosphere','Base of the crust','Depth where melting begins','location','best')

figure(7)
plot(time,vi_dim_plot*100*3600*24*365)
xlabel('Time [Gyrs]');
ylabel('Velocity [cm/yr]');

fid=fopen('CMF_sweep/Earth_melt_prod_kmGyr_CMF_0.325mmmmmmm.dat','w');
for l=1:1:length(time)
    fprintf(fid,'%8.6e %8.6e\r\n',...
        time(l),(vol_flux_plot(l)/A_planet)*1e9*3600*24*365/1e3);
end
fclose(fid);

figure(8)
plot(time,(vol_flux_plot/A_planet)*1e9*3600*24*365/1e3)
xlabel('Time [Gyrs]');
ylabel('Melt Production [km/Gyr]');

fid=fopen('CMF_sweep/Earth_heat_fluxWm2_CMF_0.325mmmmmmm.dat','w');
for l=1:1:length(time)
    fprintf(fid,'%8.6e %8.6e %8.6e\r\n',...
        time(l),ql(l),q_magma(l));
end
fclose(fid);

figure(9)
plot(time,ql,time,q_magma)
xlabel('Time [Gyrs]');
ylabel('Heat Flux [W/m^2]');
legend('Convective Heat Flux','Magmatic Heat Flux','Location','Best')

fid=fopen('CMF_sweep/Earth_heat_flowTW_CMF_0.325mmmmmmm.dat','w');
for l=1:1:length(time)
    fprintf(fid,'%8.6e %8.6e %8.6e\r\n',...
        time(l),ql(l)*A_planet/1e12,q_magma(l)*A_planet/1e12);
end
fclose(fid);

figure(10)
plot(time,ql*A_planet/1e12,time,q_magma*A_planet/1e12)
xlabel('Time [Gyrs]');
ylabel('Heat Flow [TW]');
legend('Convective Heat Flux','Magmatic Heat Flux','Location','Best')

% 
% figure(7)
% plot(time,Q_crust_tot_plot./Q_man_tot_plot)
% xlabel('Time [Gyrs]');
% ylabel('Crustal Heat production over mantle heat production (W)');
% 
% figure(8)
% plot(time,x_c_plot./x_m_plot)
% xlabel('Time [Gyrs]');
% ylabel('Crustal Heat production over mantle heat production (W/m^3)');

