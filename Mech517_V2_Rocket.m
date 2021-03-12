%Mech517_V2_Rocket
% PROBLEM 1)b. ON ASSIGNMENT 2.
% CALCULATE APOGEE OF V2 ROCKET AND GENERATE ACCELERATION, VELOCITY, AND
% ALTITUDE CURVES

%%  SET STATIC VALUES
clear 
clc
M_p = 8610; %[kg]
M_o = 12700; %[kg]
%initial weight parameters

r_o = 6.371*10^6; 
%   initial height radius of earth [m]

d = 1.6256;
%   diameter of rocket [m]
A = (pi*(d/2)^2)+(1.15824);
%   cross sectional area [m^2]

t_b = 60; %[s]
%   burn time

K = 1.4;
%   specific heat ratio air
R_air = 287;
%   ideal gas constant for air 
Isp = 250;
%   isp [s]

%%  SET INITIAL VARIABLES

g_o = 9.81;
%   initial gravity constant [m/s]
h = 0;
%   initialize height variable [m]
v = 0;
%   initialize height variable [m/s]
a = 0;
%   initialize acceleration variable [m/s^2]
M = M_o;
%   total mass of vehicle [kg]

%%  ITERATION INITIALIZE

i = 1;
t_del = .1;

while h >=0
%%  TRANSIENT PROPERTY CALCULATIONS

    g = g_o * (r_o / (r_o + h));
    %   calculate acceleration due to gravity [m/s]
   
    if h < 80,000;
        temp = 300 - .00125 * h;
    else 
        temp = 200;
    end
    
     rho=(1.2*exp((-2.9*(10^-5))*(h^1.15)));
    %   NASA earth atmosphere model for density [kg/m^3]
    c_mach = sqrt(K*R_air*(temp));
    %   speed of sound based on temp
    mach = abs(v)/c_mach;
    %   calculate mach number
    mach_stor(i) = mach;
    %   store mach number
    
    if abs(mach) <=1.1
        Cd = -1.5152*mach^4 + 3.9355*mach^3 - 2.9003*mach^2 + 0.7509*mach + 0.0953;
    elseif abs(mach) > 1.1 && abs(mach) <1.2
        Cd = .425;
    elseif abs(mach)<= 4
        Cd = -0.0294*mach^3 + 0.2735*mach^2 - 0.8528*mach + 1.09017;
    else 
        Cd = .15;
    end
    D = Cd*(1/2)*rho*v^2*A;              
    %   drag coeff and drag force from mach number [kg]
    
%%  VELOCITY CALCULATIONS

if M > (M_o-M_p)
%   during burn
v_del =(g_o * Isp * ((M_p*t_del)/(t_b*M)))-((1/M)*D*t_del)-(g*t_del);

M = M - ((M_p/t_b)*t_del);
Burn_Time = i* t_del;

elseif v > 0
%   during coast

v_del = -(1/M)*D*t_del - g*t_del;
else
%   after apogee

v_del = (1/M)*D*t_del - g*t_del;
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

%%  Integrate for acceleration curve

s = size(v_stor);
for n = 1 : (s(2) -1 )
a = (v_stor(n+1) - v_stor(n)) / t_del;
a_stor(n) = a;
end
a_stor(n+1) = a_stor(n);

%%  Plot and find maximums

Height_Reached = max(h_stor(1:4500))
V_max = max(v_stor(1:4500))
A_max = max(a_stor(1:4500))

subplot(1,3,1);
plot(t_stor(1:4850),v_stor(1:4850),'r');
xlabel('Time (t) [s]');
ylabel('Velocity [m/s]');
grid on;
title('Velocity vs Time');

subplot(1,3,2);
plot(t_stor(1:4850),h_stor(1:4850));
xlabel('Time [s]');
ylabel('Height [m]');
grid on;
title('Height vs Time');

subplot(1,3,3);
plot(t_stor(1:4850),a_stor(1:4850),'r');
xlabel('Time [s]');
ylabel('Acceleration [m/s^2]');
grid on;
title('Acceleration vs Time');