%% Bottle Rocket Problem
clear
clc
g_o = 9.81;
theta = 90 * (pi/180);
r_o = 6.371*10^6; 
%   initial height radius of earth [m]

d = 1.6256;
%   diameter of rocket [m]
A = (pi*(d/2)^2)+(1.15824);
%   cross sectional area [m^2]

Mo = 12700; %kg total mass before burn
M = Mo; %starting mass for iterations
Mp = 8610; %kg mass propellant

rho = 1.225; 

Isp = 250; %s
C = Isp*g_o;

Cd = 0.50;

K = 1.4;
%   specific heat ratio air
R_air = 287;
%   ideal gas constant for air 

V = 0;
Y_pos = 0;
X_pos = 0;
Y_v = 0;
X_v = 0;

t_b = 60; %burn time
M_dot = Mp/t_b; %kg/s flow rate

%% Iteration Initialize
i=1; %Iteration
t=0; %time set for start
t_del = .1;
M_del = M_dot * t_del; % mass change per step size

%% Burn Iteration

while Y_pos >= 0
  % Burnout at 60s
  if t >= t_b
      M_del = 0;
  end
  % Define theta based on time during flight
  
  if t < 5
      theta = 90 * (pi/180);
  elseif t < 60
      theta = 90 * (pi/180);
  else
      theta = atan(Y_v/X_v);
  end
  
  % set temperature based on height
  if Y_pos < 80,000;
        temp = 300 - .00125 * Y_pos;
  else 
        temp = 200;
  end
  
  rho=(1.2*exp((-2.9*(10^-5))*(Y_pos^1.15)));
  %   NASA earth atmosphere model for density [kg/m^3]
    
  c_mach = sqrt(K*R_air*(temp));
  %   speed of sound based on temp
  mach = abs(V)/c_mach;
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
   % Calculate drag coefficient
    
    g = g_o * (r_o / (r_o + Y_pos));
    %   calculate acceleration due to gravity [m/s]
   
    v_del = C * (M_del/M) - ((rho * A * Cd * V * V * t_del)/(2 * M)); % Velocity magnitude without gravity
    
    Y_v_del = v_del * sin(theta) - (g  * t_del); %Incorp gravity for y direction change
    X_v_del = v_del * cos(theta); % X velocity change
    
    Y_v = Y_v_del + Y_v; % current vel
    X_v = X_v_del + X_v;
    
    V = sqrt(Y_v^2 + X_v^2);
    
    Y_pos = Y_pos + (Y_v * t_del); % current position
    X_pos = X_pos + (X_v * t_del);
    
    M = M - M_del;
    
    
    
    i = i+1;
    t = i*t_del;
    t_stor(i) = t;
    
    
    mach_stor(i) = mach;
    theta_stor(i) = theta;
    V_stor(i) = V;
    Y_pos_stor(i) = Y_pos;
    X_pos_stor(i) = X_pos;
    Y_vel_stor(i) = Y_v;
    X_vel_stor(i) = X_v;
end
theta_stor(1) = 90 * (pi/180);

%% DATA PROCESSING 

s = size(V_stor);
for n = 1 : (s(2) -1 )
a = (V_stor(n+1) - V_stor(n)) / (9.81*t_del);
a_stor(n) = a;
end
a_stor(n+1) = a_stor(n);


subplot(1,2,2);
plot(t_stor,a_stor);
xlabel('Time (t) [s]');
ylabel('Acceleration [g]');
grid on;
title('Acceleration');

subplot(1,2,1);
plot(t_stor,mach_stor,'r');
xlabel('Time (t) [s]');
ylabel('Mach Number');
grid on;
title('Mach # vs Time');


subplot(2,2,2);
plot(t_stor,V_stor,'r');
xlabel('Time (t) [s]');
ylabel('Velocity [m/s]');
grid on;
title('Velocity vs Time');

subplot(2,2,1);
plot(X_pos_stor,Y_pos_stor);
xlabel('X Position [m]');
ylabel('Y Position [m]');
grid on;
title('Trajectory');

Distance = max(X_pos_stor)
Height = max(Y_pos_stor)