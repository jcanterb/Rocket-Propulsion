%% Bottle Rocket Problem
H=0;
V = 0;
g = 9.81;
Mo = 1; %kg total mass before burn
M = 1; %starting mass for iterations
Mp = 0.66666667; %kg mass propellant
A = .015; %crosssectional area
rho = 1.225; 
C = 22.09;
Cd = 0.50;
theta = 45 * (pi/180);
Y_pos = 0;
X_pos = 0;
Y_vel = 0;
X_vel = 0;
P_del = 243413;

M_dot = 8.0633; %kg/s flow rate
t_b = .08247; %burn time

%% Iteration Initialize
i=0;
t_del = .001;
M_del = M_dot * t_del; % mass change per step size

%% Burn Iteration

while M >= (Mo-Mp)
    %C = sqrt((2*P_del)/997);
   % P_del = P_del - 294.33252725;
    
    v_del = C * (M_del/M) - ((g * cos(theta)) * t_del) - ((rho * A * Cd * V * V * t_del)/(2 * M));
    
    V = v_del + V;
    Y_pos = Y_pos + (V * t_del * sin(theta));
    X_pos = X_pos + (V * t_del * cos(theta));
    Y_vel = V * sin(theta);
    X_vel = V * cos(theta);
    
    M = M - M_del;
    
    i = i+1;
    t_stor(i) = i*t_del;
    V_stor(i) = V;
    Y_pos_stor(i) = Y_pos;
    X_pos_stor(i) = X_pos;
    Y_vel_stor(i) = Y_vel;
    X_vel_stor(i) = X_vel;
end

%% Coast Iterations

while Y_pos > 0.1
    
    theta = atan(Y_vel/X_vel);
    v_del_y = - (g * t_del) - ((rho * A * Cd * V * V * t_del*sin(theta))/(2 * M));
    v_del_x = - ((rho * A * Cd * V * V * t_del*cos(theta))/(2 * M));
    
    Y_vel = v_del_y + Y_vel;
    X_vel = v_del_x + X_vel;
    Y_pos = Y_pos + (Y_vel * t_del);
    X_pos = X_pos + (X_vel * t_del);
    
    i = i+1;
    t_stor(i) = i*t_del;
    V_stor(i) = sqrt((Y_vel^2) + (X_vel^2));
    Y_pos_stor(i) = Y_pos;
    X_pos_stor(i) = X_pos;
    Y_vel_stor(i) = Y_vel;
    X_vel_stor(i) = X_vel;
end

subplot(1,2,1);
plot(t_stor,V_stor,'r');
xlabel('Time (t) [s]');
ylabel('Velocity [m/s]');
grid on;
title('Velocity vs Time');

subplot(1,2,2);
plot(X_pos_stor,Y_pos_stor);
xlabel('X Position [m]');
ylabel('Y Position [m]');
grid on;
title('Trajectory');
Distance = max(X_pos_stor)