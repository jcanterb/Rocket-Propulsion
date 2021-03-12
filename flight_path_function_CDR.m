function [Height_Reached] = flight_path_function_CDR(M_p)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%   M_p = 10.475 + 37.946 ;
%%  SET STATIC VALUES
M_l = 2.96 + 8.8 ;
M_s = 114 ; 

M_p = (M_p)/32.2;
M_f = (M_l + M_s)/32.2;
%   final rocket mass [lbm]
r_o = 3958.542 * 5280; 
%   initial height radius of earth, assume fort collins [ft]
d = 6.12/12;
%   diameter of rocket [ft]
A = (pi*(d/2)^2)+(4*.005);
%   cross sectional area [ft^2]
m_dot = 3.0375/14.594;
%   mass flow rate [slugs/s]
K = 1.4;
%   specific heat ratio air
R_air = 1716;
%   ideal gas constant for air [ft*lbf/slug*R]
C = 6204.4;
%   effective exhaust velocity [ft/s]

%%  SET INITIAL VARIABLES

g_o = 32.2;
%   initial gravity constant [ft/s]
h = 0;
%   initialize height variable [ft]
v = 0;
%   initialize height variable [ft/s]
a = 0;
%   initialize acceleration variable [ft/s^2]
M = M_f + M_p;
%   total mass of vehicle [lbm]

%%  ITERATION INITIALIZE

i = 1;
t_del = 0.005;

while h >=0
%%  TRANSIENT PROPERTY CALCULATIONS

    g = g_o * (r_o / (r_o + h));
    %   calculate acceleration due to gravity [ft/s]
    %T = 59 - (0.00356 * h);  %   [deg F]
    %P = 2116 * (((T + 459.7) / 518.6) ^ 5.256);
   % rho = P / (1718 * (T + 459.7));
    rho=(1.2*exp((-2.9*(10^-5))*(h^1.15)))*(1/515.379);
    %   NASA earth atmosphere model for density [slugs/ft^3]
    c_mach = sqrt(K*R_air*(559.6));
    %   speed of sound based on temp
    mach = abs(v)/c_mach;
    %   calculate mach number
    mach_stor(i) = mach;
    %   store mach number
    
%     if mach <=1.2
%         Cd = 0.2431*mach^3 - 0.0923*mach^2 + 0.027*mach + 0.1531;
%     elseif mach <= 5.4
%         Cd = 0.0064*mach^4 - 0.957*mach^3 + 0.5173*mach^2 -1.2229*mach + 1.278;
%     else
%         Cd = 0.15;
        
    if mach <=1.2
        Cd = 0.3472*mach^3 - 0.2381*mach^2 + 0.0188*mach + 0.1518;
    elseif mach <= 2.00
        Cd = 0.2604*mach^3 - 0.9821*mach^2 + 0.8824*mach + 0.3293;
    elseif mach <= 3.4
        Cd = 0.0063*mach^3 - 0.0318*mach^2 - 0.0086*mach + 0.3442;
    elseif mach <= 5.4
        Cd = -0.0021*mach^3 + 0.0397*mach^2 - 0.2495*mach + 0.6678;
    else
        Cd = 0.15;
    end
    D = Cd*(1/2)*rho*v^2*A;              
    %   drag coeff and drag force from mach number [lb]
    
%%  VELOCITY CALCULATIONS

if M > M_f
%   during burn
v_del =(C*(m_dot/M)*t_del)-((1/M)*D*t_del)-(g*t_del);
M = M - (m_dot*t_del);
t_b = i * t_del;
Burn_Time = t_b;
elseif v >= 0
%   during coast
v_del = -(1/M)*D*t_del - g*t_del;
else
%   after apogee
v_del = (1/M)*D*t_del-g*t_del;
end

v = v + v_del;
v_stor(i) = v;
h = h + (v * t_del);
h_stor(i) = h;
t = i * t_del;
t_stor(i) = t;
%   calc current height and velocity. store values.
i = i + 1;
end
format longG
Height_Reached = max(h_stor);
Burn_Time
s = size(v_stor);
for n = 1 : (s(2) -1 );
a = (v_stor(n+1) - v_stor(n)) / t_del;
a_stor(n) = a;
end
a_stor(n+1) = a_stor(n);

subplot(1,4,1);
plot(t_stor,v_stor,'r');
xlabel('Time (t) [s]');
ylabel('Velocity [ft/s]');
grid on;
title('Velocity vs Time');

subplot(1,4,2);
plot(t_stor,h_stor);
xlabel('Time [s]');
ylabel('Height [ft]');
grid on;
title('Height vs Time');

subplot(1,4,3);
plot(t_stor,a_stor,'r');
xlabel('Time [s]');
ylabel('Acceleration [ft/s^2]');
grid on;
title('Acceleration vs Time');

subplot(1,4,4);
plot(t_stor,mach_stor,'r');
xlabel('Time [s]');
ylabel('MACH #');
grid on;
title('MACH # vs Time');