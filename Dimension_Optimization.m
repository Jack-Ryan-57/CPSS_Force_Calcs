% Loop through brake plate angles
phi = 0:0.001:pi/2;
% Test different distances from hinge
x_a = 3:0.1:9;
% Test different distances of arm length
a = 4:0.1:8;

% Constants - area of brake plate
arc_angle_bp = 40;
L = x_a+2;
R_bt = 2.5;

% Drag data
% Last year's data for drag coefficient corresponding to max velocity
max_vel_fts = 709.0308; % Max velocity in ft/s
max_vel_mach = max_vel_fts*0.00088863; % Max mach no. 
max_vel = max_vel_fts*0.3048; % Max velocity in m/s
Cd_test = readtable('Airbrake_CD_Test.CSV');
Cd_test_vel = Cd_test{:, "Mach"};
[max_vel_tested, max_vel_index] = min(abs(max_vel_mach-Cd_test_vel)); % Find closest tested mach no. to actual max mach no.
Cd_at_max_v = Cd_test{:,"CD"}(max_vel_index); % Use maximum drag coefficient

% Find maximum angle for each dimension, return dimensions that maximize
% angle
max_angle = zeros(1, (length(x_a)*length(a)));
max_ang = 0;
index_xa = 1;
index_a = 1;
values = zeros(length(index_xa), length(index_a));
forces = zeros(length(index_xa), length(index_a));
for c = 1:length(x_a)
    for d = 1:length(a)
       max_rot = find_max_rotation(x_a(c), a(d));
       values(c, d)= (max_rot*180/pi);
       forces(c, d)= get_max_force(R_bt, x_a(c), arc_angle_bp, L(c), max_rot, Cd_at_max_v, max_vel, 4);
       if max_rot > max_ang
           max_ang = max_rot;
           index_xa = c;
           index_a = d;
       end
    end
end
fprintf("Maximum rotation angle: %d degrees\n", max_ang*180/pi);
fprintf("Preferred mounting distance from hinge: %d inches\n", x_a(index_xa));
fprintf("Preferred arm length: %d inches\n", a(index_a));
fprintf("Maximum force on motor with given dimensions: %d lbf", forces(index_xa, index_a));

function [max_phi] = find_max_rotation(x_a, a)
    R_bt = 2.5;
    phi = 0:0.01:pi/2;
    r_ca = sqrt(x_a.^2+R_bt^2-2.*x_a.*R_bt.*cos(pi/2+phi));
    theta = asin((R_bt+x_a.*sin(phi))./(a));
    [d, ~] = solve_quadratic(1, -2*a.*cos(theta), (a.^2-r_ca.^2)); % Solve for d
    [~, ind] = find_complex(d);
    max_phi = phi(ind);  % Maximum rotation for given x_a and a
end

% Quadratic formula
function [pos_root, neg_root] = solve_quadratic(a, b, c)
    pos_root = (-b+sqrt(b.^2-4.*a.*c))/(2.*a);
    neg_root = (-b-sqrt(b.^2-4.*a.*c))/(2.*a);
end

% Find the first complex value in array of values
function [max, ind] = find_complex(arr)
    ind = 1;
    for c = 1:length(arr)
        if ~isreal(arr(c))
            ind = c-1;
            break;
        end
    end 
    max = arr(ind);
end

function force = get_max_force(R_bt, x_a, arc_angle_bp, L, rotation, Cd_at_max_v, max_vel, num_brakes)
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
    force = num_brakes.*F_on_actuator;
end
