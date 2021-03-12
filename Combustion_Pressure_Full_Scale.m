%% Tank Geometry

f_tank_id = 2.62; % in
o_tank_id = 5.62; % in
f_tank_thick = .13; % in
tank_height = 55; % in originally 91.88
init_o_headspace = 3.5; % in
init_o_level = tank_height - init_o_headspace; % in
init_f_level = init_o_level - 2.7; % in

%% Tank Initial Properties

init_tank_temp_F = 64.3; % deg F
init_tank_temp_C = (init_tank_temp_F - 32)*(5/9); % deg C

f_MW = 46;
f_dens = 789;
f_dens_lb = f_dens /16.018; % lbm/ft^3

o_MW = 44;
o_dens = -0.1879*(init_tank_temp_C^2) - 1.8974*init_tank_temp_C + 899.47; % kg/m3 From polynomial for liquid density
o_dens_lb = o_dens/16.018; % lbm/ft^3

f_vol = init_f_level * pi * 0.25 * (f_tank_id^2); % in^3
o_vol = (init_o_level * pi * 0.25 * (o_tank_id^2)) - (init_o_level * pi * 0.25 * ((f_tank_id +(f_tank_thick * 2))^2)); % in^3

f_mass_lb = (f_vol/1728) * f_dens_lb; % lb
o_mass_lb = (o_vol/1728) * o_dens_lb; % lb

del_P_fpiston_lb = 25; % psi
del_P_fpiston = del_P_fpiston_lb *(101325/14.7); % Pa

tank_P = 0.0098 * (init_tank_temp_C^2) + (0.7715*init_tank_temp_C) + 31.265; % bar
tank_P_lb = tank_P * 14.5038; %psi

%% Injector Specs

inject_D = 0.08; %in
inject_area = pi * .25 * (inject_D^2); %in^2
inject_coeff = .45; 
n_inject = 6;
inject_area_metric = inject_area*(2.54/100)^2; %m^2

annulus_OD = .76;  %in
annulus_ID = .45;  %in
annulus_area = pi * ((annulus_OD/2)^2-(annulus_ID/2)^2); %in^2
annulus_coeff = .175;
n_annulus = 1;
annulus_area_metric = annulus_area*(2.54/100)^2; %m^2

%% Nozzle Specs

dt = 3; %in
At = pi*.25*(dt^2); %in^2
At_metric = At*((2.54/100)^2); %m^2
Pe_Pc = 14/350;
Pa = 850* (14.7/1013); % psi
gamma_approx = 1.2;
G =(gamma_approx^0.5)*(2/(1+gamma_approx))^((gamma_approx+1)/2/(gamma_approx-1));
Ae_At = G/((Pe_Pc^(1/gamma_approx))*(2*gamma_approx*(1-Pe_Pc^((gamma_approx-1)/gamma_approx))/(gamma_approx-1))^0.5);
Ae = Ae_At * At;
de = 4.5;

C_star_eff = 1;
Cf_eff = .85;
%% Iterate to find pressure static
es = .00001;  % set error limit
ea= 100;    % set error at 100% to start
P_c = 300; %psi
i = 0;  %count iterations

while ea >= es
    del_P = (tank_P_lb - P_c)*(101325/14.7); %Pa
    m_dot_f = inject_coeff * n_inject * inject_area_metric * sqrt(2 * f_dens * (del_P - del_P_fpiston)); %kg/s
    m_dot_o = annulus_coeff * n_annulus * annulus_area_metric * sqrt(2 * o_dens * del_P); %kg/s

    m_dot_tot = m_dot_o + m_dot_f; % kg/s

    o_f = m_dot_o / m_dot_f; 

    C_star = ((12.456 * (o_f^3)) - (234.37 * (o_f^2)) + (1352.2 * o_f) + 2780)* C_star_eff; %ft/s
    
    P_old = P_c; %psi
    
    P_c = (C_star*0.3048*(m_dot_tot)/At_metric)*14.7/101325; %psi

    ea = abs(((P_old - P_c)/P_c)*100); % percent error

    i=i+1;
end
gamma = (0.00023*o_f^5)-(0.006436*o_f^4)+(0.06815*o_f^3)-(0.33097*o_f^2)+(0.6821*o_f)+0.78795;

C_f_o = sqrt((2*gamma^2/(gamma-1)*(2/(gamma+1))^((gamma+1)/(gamma-1))*(1-(Pe_Pc^((gamma-1)/gamma)))));
C_f =(C_f_o+(Pe_Pc-(Pa/P_c))*Ae_At)*Cf_eff; %ft/s

F = C_f * P_c * At; %lbf
Isp = (C_f * C_star)/32.2; % s
C = C_f * C_star;

%%Burnout and Total Impulse ESTIMATE

o_burnout_time = o_mass_lb /(m_dot_o*2.20462);
f_burnout_time = f_mass_lb / (m_dot_f*2.20462);

total_impulse = o_burnout_time * F;

%% Display results

format shortG
disp('     Chamber P      Tank P       Mdot F        Mdot O     Mdot Total')
disp('       [psi]         [psi]       [kg/s]        [kg/s]       [kg/s]')
x =[P_c tank_P_lb m_dot_f m_dot_o m_dot_tot];
disp(x)

disp('        O/F      Eff Velocity      F           Isp        O_burnout    F_burnout    Tot Impulse')
disp('        [-]         [ft/s]       [lbf]         [s]           [s]          [s]        [lbf*s]')
x =[o_f C F Isp o_burnout_time f_burnout_time total_impulse];
disp(x)