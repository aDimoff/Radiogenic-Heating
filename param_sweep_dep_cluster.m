close all
clear variables

% Dorn uses 1000 ppm in the mantle, which is 1e-3 wt %; so same as my low
% end case. They then have only 10% of this melt erupt at surface, but also
% don't deplete the mantle of CO2
%
% Also, Dorn only outgasses 10% of melt; I could look into this for the LHS
% 3844 paper, and do the same; rest of CO2 trapped in crust. Would bring
% back a metamorphic outgassing flux; but no real recycling, as CO2 would
% still go to the atmosphere and stay there

% Single model, includes delamination and crust decarbonation
% New melt heat loss treatment
% Includes partition coefficient for CO2 degassing

% Mass of the Earth
M_E=5.97e24; % in kg 

% Set the planet mass (in Earth masses) and code calculates planet radius,
% mantle radius, density, etc based on structure models from Valencia et al
% 2007
mass=1.0; 
Rp=6378100.0*mass^0.27;
rho=4450.0*mass^0.2;
g=mass*M_E*6.67e-11/Rp^2;
d=2890000*mass^0.28;
Rc=Rp-d; % Core radius

Ts=273.0; % Surface temperature in K

% Material parameters
Ev=300000.0; % Activation energy for viscosity, J/mol
Rg=8.314; % Gas constant
alpha=3e-5; % Themal expansion coefficient, 1/K
kappa=1e-6; % Thermal diffusivity m^2/s
rho_m=3300.0; % Density at surface pressure kg/m^3
dphi_dP=0.15e-9; % Chang in melt fraction with pressure 1/Pa
Cp=1250; % Heat capacity, J kg^-1 K^-1
k=5; % Thermal conductivity, W m^-1 K^-1
C1=0.5; % heat flow scaling law constant
L_m=600*1e3; % latent heat J/kg
dT_dP=2e-8; % adiabatic gradient in mantle

% This parameter sets the fraction of devolatilized carbon that is trapped 
% and still recycles into mantle. Used for testing in initial papers but
% since then always set to 0
f=0;

V_man=(4/3)*pi*(Rp^3-Rc^3); % Mantle volume for the modeled planet
V_man_Earth=(4/3)*pi*(6378100^3-3488100^3); % Earth's mantle volume
m_man=V_man*rho; % Mass of the mantle 
m_man_Earth=3.9e24; % mass of Earth's mantle in kg

% Intrinsic heat production for each of the 4 key HPEs
Q_U238=9.17e-5; % W/kg of U238
Q_U235=5.75e-4; % W/kg of U235
Q_Th=2.56e-5; % W/kg of Th
Q_K=2.97e-5; % W/kg of K

% Half-lives
tau_half_U238=4.46e9; % years
lambda_U238=-tau_half_U238/log(1/2);
tau_half_U235=7.04e8; % years
lambda_U235=-tau_half_U235/log(1/2);
tau_half_Th=1.4e10; % years
lambda_Th=-tau_half_Th/log(1/2);
tau_half_K=1.26e9; % years
lambda_K=-tau_half_K/log(1/2);

% Determine starting radiogenic heating rates for each radionuclide
% This extrapolates present day bulk silicate Earth values back to the
% beginning of Earth history 
%
% For stars with different ages you would want to extrapolate back based on
% the age of the star
Earth_U238=0.022*0.9927*exp(4.5e9/lambda_U238); %ppm
Earth_U235=0.022*0.0072*exp(4.5e9/lambda_U235); %ppm
Earth_Th=0.083*exp(4.5e9/lambda_Th); %ppm
Earth_K=261*(0.0117/100)*exp(4.5e9/lambda_K); %ppm

% Lognormal distribution for the HPEs based on Hypatia catalog
% Unterborn et al, 2022
mu_U238=1.0;
sigma_U238=0.25;
mu_U235=1.0;
sigma_U235=0.25;
mu_Th=1.0;
sigma_Th=0.3;
mu_K=1.0;
sigma_K=0.3;

% End time to run the model to (in seconds)
t_end=50e9*3600*24*365;

Ae=4*pi*Rp^2; % surface area of planet
Ae0=4*pi*6378100^2; % surface area of Earth
m_bar_co2=44/1000; % molar mass of CO2
m_bar_w=18/1000; % molar mass of H20
plume_vol_flux=0;
rho_r=2800;% Used as melt density, NOT basalt density in this code
Fdstar=6e12/(3600*24*365); % Present day CO2 degassing flux on Earth (mol/s)

% These parameters are to increase the solidus temperature due to chemical
% depletion of the mantle
dcr_crit=0.2*V_man/Ae;
DT_sol_dep=150; 

% Supply limited weathering rate equation
% First need moles of CaO, FeO, MgO, K2O, Na2O per 100 g 
% Weight percents: MgO 7.58; CaO 11.39; FeO 10.43; Na2O 2.79; K2O 0.16
% Weight percents from Gale et al, 2013
% Mol/100 g: MgO 0.19; CaO 0.2; FeO 0.145; Na2O 0.045; K2O 0.0017
weath_demand=10*(0.19+0.2+0.145+0.045+0.0017); % mol/kg of CO2 drawdown

% These are scaling factors to model galactic chemical evolution in
% Unterborn et al 2022
% Setting them to 1 means using present day values, i.e. not considering
% galactic chemical evolution
ga_scale_K=1.0;
ga_scale_Th=1.0;
ga_scale_U235=1.0;
ga_scale_U238=1.0;

ex_time=1;
rng('shuffle');
fid=fopen('Earth_dep_degassing_lifetimes','w');

% This is the starting CO2 budget of the planet in mol
% Can set the value directly or set a value that scales with mantle mass 

%C_tot=5e21*(m_man/m_man_Earth); 
C_tot=5e21; 

for j=1:50000
    % Sample uniform distributions of ref viscosity, C budget, initial T
    mun=10^unifrnd(9.3444,13.3444); % sample from mun=4e10 to 4e12;
    %C_tot=10^unifrnd(18.5,23.5)*(m_man/m_man_Earth);
    Ti_init=unifrnd(1700,2000);
    muref=mun*exp(Ev/(Rg*1623));
    
    % Sample from log normal distributions of HPEs
    U238=max(0,normrnd(mu_U238,sigma_U238));
    U235=max(0,normrnd(mu_U235,sigma_U235));
    Th=max(0,normrnd(mu_Th,sigma_Th));
    K=max(0,normrnd(mu_K,sigma_K));    
    
    Q0_U238=ga_scale_U238*m_man*Q_U238*U238*Earth_U238*(1/1e6);
    Q0_U235=ga_scale_U235*m_man*Q_U235*U235*Earth_U235*(1/1e6);
    Q0_Th=ga_scale_Th*m_man*Q_Th*Th*Earth_Th*(1/1e6);
    Q0_K=ga_scale_K*0.5*m_man*Q_K*K*Earth_K*(1/1e6);
    Q0_tot=Q0_U238+Q0_U235+Q0_Th+Q0_K;
    
    Q0_U238_save(j)=Q0_U238;
    Q0_U235_save(j)=Q0_U235;
    Q0_Th_save(j)=Q0_Th;
    Q0_K_save(j)=Q0_K;
    Q0_tot_save(j)=Q0_tot;
    
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
    % [T,Y]=ode15s(@(t,y) thermal_evol_delam_melt_full_Qpartition_atm2(t,y,L_m,1,mun,h_phi),...
    %     [0 t_end],[Ti_init 0 0 Rman_init 0 0 Q0_forward delta_init Raoc_init],options);
    [T,Y,Te,Ye,ie]=ode15s(@(t,y) thermal_evol_delam_melt_full_Qpartition_atm2_dep(t,y,L_m,1,mun,Ts,Rp,rho,g,...
        d,Rc,dT_dP,dcr_crit,DT_sol_dep),[0 t_end],[Ti_init 0 0 1.0*C_tot 0 0.0*C_tot Q0_U238 ...
        delta_init 0 Q0_U235 0 Q0_Th 0 Q0_K],options);
    toc
    time(1:length(T))=T(1:length(T))/(1e9*3600*24*365);
    Ti(1:length(T))=Y(1:length(T),1);
    crust(1:length(T))=Y(1:length(T),2);
    crust2(1:length(T))=Y(1:length(T),5);
    Rman(1:length(T))=Y(1:length(T),4);
    Rcrust(1:length(T))=Y(1:length(T),6);
    delta(1:length(T))=Y(1:length(T),8);
    u238_crust(1:length(T))=Y(1:length(T),3);
    u235_crust(1:length(T))=Y(1:length(T),9);
    th_crust(1:length(T))=Y(1:length(T),11);
    K_crust(1:length(T))=Y(1:length(T),13);
    u238_man(1:length(T))=Y(1:length(T),7);
    u235_man(1:length(T))=Y(1:length(T),10);
    th_man(1:length(T))=Y(1:length(T),12);
    K_man(1:length(T))=Y(1:length(T),14);
    
    % Sorting events
    supply_lim_start=10;
    fdegas1=0;
    fdegas2=0;
    fdegas_end=0;
    fdegas_start=0;
    for i=1:1:length(Te)
        if ie(i) == 1
            fdegas1=Te(i)/(3600*24*365*1e9);
        elseif ie(i) == 2
            fdegas2=Te(i)/(3600*24*365*1e9);
        elseif ie(i) == 3
            fdegas_end=Te(i)/(3600*24*365*1e9);
        elseif ie(i) == 4
            supply_lim_start=Te(i)/(3600*24*365*1e9);
        elseif ie(i) == 5
            fdegas_start=Te(i)/(3600*24*365*1e9);
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
        %d_melt_plot(l)=min(d_melt,d_dense);
        if (P_melt > P_lid)
            phi=0.5*(P_melt1-P_lid)*dphi_dP;
        else
            phi=0;
        end
        vi_dim=(kappa/d)*0.05*(Rai/theta)^(2/3);
        %     vi_dim=(kappa/d)*0.05*(Rai*min(1/theta,...
        %         C1^(-3/4)*theta^(1/3)*(theta^(1/3)*(theta-1))^(-1/4)*(Y(2)/(d*Rai))^(1/4)*(d/Y(2))))^(2/3);
        vol_flux=17.8*3.14*Rp^2.0*vi_dim*phi*(min(d_melt,d_dense)-delta(l))*(1/d);
        %     vol_flux2=64*3.14*Rp^2.0*vi_dim*(d_melt-delta)*(1/d); % Upwelling mantle into melt zone flux
        
        % For partitioning of CO2 into melt
        frac2=(1-(1-phi)^(1/1e-4));
        
        % Temp at base of the lid
        Tl=Ti(l)-2.5*(Rg*Ti(l)^2)/Ev;
        
        Q_crust_tot=u238_crust(l)+u235_crust(l)+th_crust(l)+K_crust(l);
        Q_man_tot=u238_man(l)+u235_man(l)+th_man(l)+K_man(l);
        % Concentration of heat producing elements (watts per meter cubed) in crust
        % and mantle
        if (crust2(l)==0)
            x_c=0;
        else
            x_c=Q_crust_tot/crust2(l);
        end
        x_m=Q_man_tot/V_man;
        
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
    if max(Fdegas+F_dcarb-Fws1)*3600*24*365>1e14
        supply_lim2(j)=1;
    else
        supply_lim2(j)=0;
    end

    degas_sum=0;
    i_degas=find(time<fdegas_end & time>fdegas_start);
    if isempty(i_degas)
        degas_avg=0;
    else
        for l=1:length(i_degas)-1
            degas_sum=degas_sum+(1/2)*(Fdegas(i_degas(l))+F_dcarb(i_degas(l))+...
                Fdegas(i_degas(l+1))+F_dcarb(i_degas(l+1)))*(time(i_degas(l+1))-time(i_degas(l)))*(3600*24*365*1e9);
        end
        degas_avg=degas_sum/(time(i_degas(end))*(3600*24*365*1e9)-time(i_degas(1))*(3600*24*365*1e9));
        degas_avg=degas_avg*3600*24*365;
    end
    
    fprintf(fid,'%8.6e %8.6e %8.6e %8.6e %8.6e %8.6e %8.6e %8.6e %8.6e %8.6e %8.6e %8.6e %8.6e %8.6e %8.6e %8.6e %8.6e\n',...
        mun,muref,C_tot,Ti_init,delta_init,Q0_tot_save(j),Q0_U238_save(j),Q0_U235_save(j),Q0_Th_save(j),...
        Q0_K_save(j),fdegas1,fdegas2,fdegas_start,fdegas_end,supply_lim_start,supply_lim2(j),degas_avg);
    
    clear time Ti crust crust2 Rman Rcrust delta u238_crust u235_crust th_crust K_crust u238_man u235_man th_man K_man
    clear Fws1 Fdegas F_dcarb i_degas
    clear T Y Te Ye ie 
    
    j
end




