% VORTICITY/STREAM FUNCTION LID DRIVEN CAVITY FLOW SOLVER BY JOE MOLVAR
clear; close all;
% Given Parameters
Nx = 52;
L = 1;
Wall_Velocity = 2.75; % Nodes X; Domain Size; Wall Velocity
rho = 1;
mu = 0.01; % Density; Dynamic Viscosity;
dt = 0.001;
maxIt = 50000;
maxe = 1e-7; % Time Step; Max iterations; Max error
% 1D Grid Setup
Ny = Nx;
h = L / (Nx - 1);
x = 0:h:L;
y = 0:h:L;
im = 1:Nx-2; i = 2:Nx-1; ip = 3:Nx; jm = 1:Ny-2; j = 2:Ny-1; jp = 3:Ny;
% Initialize Matrices
Vo = zeros(Nx, Ny);
St = zeros(Nx, Ny);
Vop = zeros(Nx, Ny);
u = zeros(Nx, Ny);
v = zeros(Nx, Ny);
u(2:Nx-1, Ny) = Wall_Velocity;
Vop = Vo;
% Solve Loop (Similar to Gauss-Siedel Method)
for iter = 1:maxIt
    % Boundary Conditions
    Vo(:, Ny) = -2 * St(:, Ny-1) / (h^2) - Wall_Velocity * 2 / h; % Top
    Vo(:, 1) = -2 * St(:, 2) / (h^2); % Bottom
    Vo(1, :) = -2 * St(2, :) / (h^2); % Left
    Vo(Nx, :) = -2 * St(Nx-1, :) / (h^2); % Right
    % Partially Solve Vorticity Transport Equation
    s_x1 = max(u(i, j), 0);
    s_x2 = max(-u(i, j), 0);
    s_y1 = max(v(i, j), 0);
    s_y2 = max(-v(i, j), 0);
    % SOR Calculation
    OMEGA = 1.5;
    Vo(i, j) =Vop(i, j) + (OMEGA ./ (1 - 4*dt*(mu/(rho*h^2) - dt*(s_x1 ./ h) -...
                dt*(s_x2 ./ h) - dt*(s_y1 ./ h) - dt*(s_y2 ./ h))) .* ((mu/(rho*h^2)) * (Vop(ip, j) + Vop(im, j) + Vop(i, jp) + Vop(i, jm) - 4 * Vop(i, j)) - ...
        (s_x1 ./ h) .* (Vop(i, j) - Vop(im, j)) + (s_x2 ./ h ).* (Vop(ip, j) - Vop(i, j)) - ... 
         (s_y1 ./ h) .* (Vop(i, j) - Vop(i, jm)) + (s_y2 ./ h) .* (Vop(i, jp) - Vop(i, j))) * dt);
    % Partially Solve Elliptical Vorticity Equation for Stream Function
    St(i, j) = (Vo(i, j) * h^2 + St(ip, j) + St(i, jp) + St(i, jm) + St(im, j)) / 4;
    u(i, j) = (St(i, jp) - St(i, jm)) / (2 * h);
    v(i, j) = (-St(ip, j) + St(im, j)) / (2 * h);
    % Check for Convergence
    if iter > 10
        error = max(max(abs(Vo - Vop)));
        if error < maxe
            break;
        end
    end
    Vop = Vo;
end
% Plots
cm = hsv(ceil(100 / 0.7));
cm = flipud(cm(1:100, :));
figure(1);
contourf(x, y, u', 23, 'LineColor', 'none');
title('U-velocity');
xlabel('x-location');
ylabel('y-location');
% contourf(x, y, v', 23, 'LineColor', 'none');
% title('V-velocity');
% xlabel('x-location');
% ylabel('y-location');
axis equal;
axis([0 L 0 L]);
colormap(cm);
colorbar('westoutside');
figure(2);
plot(y, u(round(Ny/2), :));
title('Centerline x-direction velocity');
xlabel('y/L');
ylabel('u/U');
axis square;
xlim([0 L]);
grid on;
N = 1000;
xstart = L * rand(N, 1);
ystart = L * rand(N, 1);
[X, Y] = meshgrid(x, y);
figure(3);
h = streamline(X, Y, u', v', xstart, ystart, [0.1, 200]);
title('Stream Function');
xlabel('x-location');
ylabel('y-location');
axis equal;
set(h, 'color', 'k');
% a and b
%Find indices for the middle vertical and horizontal lines

%Define the midpoints for vertical and horizontal midlines
midpoint_x = round(Nx / 2);
midpoint_y = round(Ny / 2);

% Extract the vertical velocity profile (v) along the midlines
v_vertical_midline = v(midpoint_x, :);
v_horizontal_midline = v(:, midpoint_y);

% Plot the vertical velocity profile along the vertical midline
figure;
plot(y, v_vertical_midline, 'b-', 'LineWidth', 2);
xlabel('y');
ylabel('v');
title('Vertical Velocity Profile along Vertical Midline');
grid on;
ylim([min(v(:)), max(v(:))]);

% Plot the vertical velocity profile along the horizontal midline
figure;
plot(x, v_horizontal_midline, 'r-', 'LineWidth', 2);
xlabel('x');
ylabel('v');
title('Vertical Velocity Profile along Horizontal Midline');
grid on;
ylim([min(v(:)), max(v(:))]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Extract the vertical velocity profile (v) along the midlines
u_vertical_midline = u(midpoint_x, :);
u_horizontal_midline = u(:, midpoint_y);

% Plot the vertical velocity profile along the vertical midline
figure;
plot(y, u_vertical_midline, 'b-', 'LineWidth', 2);
xlabel('y');
ylabel('u');
title('horizental Velocity Profile along Vertical Midline');
grid on;
ylim([min(u(:)), max(u(:))]);

% Plot the vertical velocity profile along the horizontal midline
figure;
plot(x, u_horizontal_midline, 'r-', 'LineWidth', 2);
xlabel('x');
ylabel('u');
title('horizental Velocity Profile along Horizontal Midline');
grid on;
ylim([min(u(:)), max(u(:))]);
%**********************part c
 stress_top(2,:) = mu*((u(i,Ny)-u(i,Ny-1))/h + (v(ip,Ny)-v(i,Ny))/h)
 stress_bottom(2,:) = mu*((u(i,2)-u(i,1))/h + (v(ip,1)-v(i,1))/h)

%**********************part c
%  stress_top(2,:) = mu*((u(i,Ny)-u(i,Ny-1))/h + (v(ip,Ny)-v(i,Ny))/h)
%  stress_bottom(2,:) = mu*((u(i,2)-u(i,1))/h + (v(ip,1)-v(i,1))/h)
%  %%% PLOT STRESS
% figure;
% subplot(2,1,1);
% plot(x(2:end-1), stress_top);
% title('Stress on the Top Wall');
% xlabel('x');
% ylabel('Stress');
% 
% subplot(2,1,2);
% plot(x(2:end-1), stress_bottom);
% title('Stress on the Bottom Wall');
% xlabel('x');
% ylabel('Stress');

