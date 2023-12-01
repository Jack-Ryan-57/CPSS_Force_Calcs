% CPSS Avionics Team 2023 
% Actuator minimum force rating calculations
% Jack Ryan

% Goal:
% Find maximum drag force on actuator to predict minimum necessary force 
% rating for actuator-driven airbrake system

% Constants
%---------------------
% Test values - subject to change
%---------------------
R_bt = 2.5; % Radius of body tube
max_angle = 79; % Maximum angle of air brake deployment
arc_angle_bp = 40; % Arc angle of brake plate (theta) in degrees
L = 7.6; % Length of brake plate in inches
len_arm = 8; % length of arm connecting actuator to brake plate
x_a = 5.6; % Distance from brake plate hinge to intersection with arm (in inches)
num_brakes = 4; % Number of brake plates on body tube
screw_pitch = 0.25;


% Last year's maximum predicted velocity 
max_vel_fts = 709.0308; % Max velocity in ft/s
max_vel_mach = max_vel_fts*0.00088863; % Max mach no. 
max_vel = max_vel_fts*0.3048; % Max velocity in m/s
% Last year's data for drag coefficient corresponding to max velocity
Cd_test = readtable('Airbrake_CD_Test.CSV');
Cd_test_vel = Cd_test{:, "Mach"};
[max_vel_tested, max_vel_index] = min(abs(max_vel_mach-Cd_test_vel)); % Find closest tested mach no. to actual max mach no.
Cd_at_max_v = Cd_test{:,"CD"}(max_vel_index); % Use maximum drag coefficient

rotation = deg2rad(0:0.1:max_angle); % Brake plate rotation (phi) in radians

% Brake plate rotation to linear actuator displacement
[d_T, ~] = solve_quadratic(1, -2*len_arm*cos(asin(R_bt/len_arm)), (len_arm^2-x_a^2-R_bt^2)); % Distance from retracted actuator to center of body tube at point of hinge (in inches)
r_ca = sqrt(x_a^2+R_bt^2-2*x_a*R_bt.*cos(pi/2+rotation)); 
arm_angle = asin((R_bt+x_a.*sin(rotation))./len_arm); % Angle that arm makes with vertical axis
[d, ~] = solve_quadratic(1, -2*len_arm.*cos(arm_angle), (len_arm.^2-r_ca.^2)); % Solve for d
delta_x = d_T - d; % Linear displacement angle as a function of brake plate rotation

% Exposed area vs. rotation
w = 2*R_bt*sind(arc_angle_bp/2);
A_exposed_in = L*w.*sin(rotation); % Exposed surface area in in^2
A_exposed = A_exposed_in.*0.00064516; % Exposed surface area in m^2

% Drag force vs. Exposed Surface Area
air_dens = 1.225; % Density of air at sea level (theoretical max)
Fd_max = 0.224809*0.5*air_dens*Cd_at_max_v*max_vel^2.*(A_exposed); % Maximum drag force (in lbf)

% Drag force to force on arm
x_drag = 0.5*L.*sin(rotation); % Horizonal location of drag on brake
F_on_actuator = x_drag.*(Fd_max)./x_a;

% Downward force to minimum motor torque
torque_req = num_brakes*max(F_on_actuator)*screw_pitch/(4*pi)*1.15212;

% Plot figures for drag force on airbrake
figure(1);
subplot(2, 1, 1);
plot(delta_x, Fd_max);
title("Actuator displacement vs. Force of Drag");
xlabel("Actuator displacement (in)");
ylabel("Drag force on airbrake (lbf)");
subplot(2, 1, 2);
plot(rad2deg(rotation), Fd_max);
title("Rotation of brake plate vs. Force of Drag");
xlabel("Rotation (degrees)");
ylabel("Drag force on airbrake (lbf)");

% Plot figures for force on actuator
figure(2);
subplot(2, 1, 1);
plot(delta_x, F_on_actuator);
title("Linear displacement of actuator vs. drag force on actuator");
xlabel("Linear displacement (in)");
ylabel("Force on actuator (lbf)");
subplot(2, 1, 2);
plot(rad2deg(rotation), F_on_actuator);
title("Brake plate rotation vs. drag force on actuator");
xlabel("Rotation (degrees)");
ylabel("Force on actuator (lbf)");

% Print values of interest
fprintf('Maximum force from single plate: %d lbs\n', max(F_on_actuator));
fprintf('Maximum force on actuator: %d lbs\n', num_brakes*max(F_on_actuator));
fprintf('Maximum exposed surface area: %d in^2\n', num_brakes*max(A_exposed_in));
fprintf('Expected minimum torque: %d kg*cm\n', torque_req);

% Quadratic formula
function [pos_root, neg_root] = solve_quadratic(a, b, c)
    pos_root = (-b+sqrt(b.^2-4.*a.*c))/(2.*a);
    neg_root = (-b-sqrt(b.^2-4.*a.*c))/(2.*a);
end
