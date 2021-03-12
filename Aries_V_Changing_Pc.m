%%  SET STATIC VALUES
clear 
clc

M_p = 24/2.205; %[kg]
M_o = 134/2.205; %[kg]
%initial weight parameters

r_o = (6.371*10^6)+1400.55; 
%   initial height radius of earth [m]

d = 6/39.37;
%   diameter of rocket [m]
A = (pi*(d/2)^2);
%   cross sectional area [m^2]

A_t = .001781;
% throat area ( math in EXCEL) [m^2]

Ae_At = 3;
% Area ratio, fixed

C_star = 1601;
% Cstar [m/s]

C_fo = 1.351229;
% nonchanging Cf piece

Pe_Pc = 0.070183;
%Pressure Ratio

K = 1.16;
%   specific heat ratio air
R_air = 287;
%   ideal gas constant for air 


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
t_del = .01;

while h >=0
%%  TRANSIENT PROPERTY CALCULATIONS

    g = g_o * ((r_o / (r_o + h))^2);
    %   calculate acceleration due to gravity [m/s]
   
    if h < 80,000;
        temp = 300 - .00125 * h;
    else 
        temp = 200;
    end
    
     rho=(1.2*exp((-2.9*(10^-5))*(h^1.15)));
    %   NASA earth atmosphere model for density [kg/m^3]
    c_mach = sqrt(1.4*R_air*(temp));
    %   speed of sound based on temp
    mach = abs(v)/c_mach;
    %   calculate mach number

    
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
%%  MOTOR PARAMETER CALCULATIONS
t = i * t_del;

Pc = (0.0725*t^5 - 1.5035*t^4 + 11.665*t^3 - 41.482*t^2 + 56.108*t + 183.54) * 6895 ;% calculate Pc as function of time and convert to Pascals

Pa = 101325 * ((1 - (0.0000225577 * h))^5.25588); % pa equation given in problem statement

m_dot = (Pc * A_t)/C_star;

m_del = m_dot * t_del;

C_f = C_fo + ((Pe_Pc - (Pa/Pc)) * Ae_At);

C = C_star * C_f;

%%  VELOCITY CALCULATIONS

if M > (M_o-M_p)
%   during burn
v_del =(C * (m_del/M)) - ((1/M)*D*t_del)-(g*t_del);

M = M - m_del;

BURNTIME = i*t_del;
tb_stor(i) = BURNTIME;
m_dot_stor(i) = m_dot;
Isp_stor(i) = C/9.81;                               % STORING BURN SPECIFIC PARAMETERS
Cf_stor(i) = C_f;
Pe_stor(i) = Pe_Pc * Pc;
Thrust_stor(i) = m_dot * C;


elseif v > 0
%   during coast

v_del = -(1/M)*D*t_del - g*t_del;
else
%   after apogee

v_del = (1/M)*D*t_del - g*t_del;
end

%%CALCULATE VARIOUS PARAMETERS FOR PRINTING


Pa_stor(i) = Pa;
mach_stor(i) = mach;


v = v + v_del;
v_stor(i) = v;
h = h + (v * t_del);
h_stor(i) = h;

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

%%  PLOT FIGURE 1

figure(1)
subplot(1,4,2);
plot(t_stor,v_stor,'r');
xlabel('Time [s]');
ylabel('Velocity [m/s]');
grid on;
title('Velocity vs Time');

figure(1)
subplot(1,4,1);
plot(t_stor,h_stor,'r');
xlabel('Time [s]');
ylabel('Height [m]');
grid on;
title('Height vs Time');

figure(1)
subplot(1,4,4);
plot(t_stor,a_stor,'r');
xlabel('Time [s]');
ylabel('Acceleration [m/s^2]');
grid on;
title('Acceleration vs Time');

figure(1)
subplot(1,4,3);
plot(t_stor,mach_stor,'r');
xlabel('Time [s]');
ylabel('Mach Number');
grid on;
title('Mach # vs Time');

%%  PLOT FIGURE 2

figure(2)
subplot(1,3,1);
plot(t_stor,Pa_stor,'r');
xlabel('Time [s]');
ylabel('Pa [Pa]');
grid on;
title('Atm Pressure vs Time');

figure(2)
subplot(1,3,2);
plot(tb_stor,Pe_stor,'r');
xlabel('Time [s]');
ylabel('Pe [Pa]');
grid on;
title('Exhaust Pressure vs Time');

figure(2)
subplot(1,3,3);
plot(tb_stor,Cf_stor,'r');
xlabel('Time [s]');
ylabel('C_f');
grid on;
title('Thrust Coeff vs Time');

%%  PLOT FIGURE 3

figure(3)
subplot(1,3,1);
plot(tb_stor,m_dot_stor,'r');
xlabel('Time [s]');
ylabel('Mass Flow Rate [kg/s]');
grid on;
title('Mass Flow Rate vs Time');

figure(3)
subplot(1,3,2);
plot(tb_stor,Thrust_stor);
xlabel('Time [s]');
ylabel('Thrust [N]');
grid on;
title('Thrust vs Time');

figure(3)
subplot(1,3,3);
plot(tb_stor,Isp_stor,'r');
xlabel('Time [s]');
ylabel('Isp [s]');
grid on;
title('Isp vs Time');