clear; close all
%%**************steady state equation*******************
%%% GIVENS
Nx = 52; L = 1; Wall_Velocity = 10; % Nodes X; Domain Size; Velocity
rho = 1; mu = 0.01;  k=0.3; % Density; Dynamic Viscosity;
 maxIt = 50000; maxe = 1e-7; % Time Step; Max iter; Max error
%%% SETUP 1D GRID
Ny = Nx; h = L / (Nx - 1); x = 0:h:L; y = 0:h:L;
im = 1:Nx-2; i = 2:Nx-1; ip = 3:Nx; jm = 1:Ny-2; j = 2:Ny-1; jp = 3:Ny;
%%% PRELOCATE MATRICES
Vo = zeros(Nx, Ny); St = Vo; Vop = Vo; u = Vo; v = Vo; T = 300 * ones(Nx, Ny); % Initial temperature

%%% VELOCITY ON THE UPPER WALL (NO-SLIP CONDITION)
u(2:Nx-1, Ny) = Wall_Velocity;

%%% SOLVE LOOP SIMILAR TO GAUSS-SEIDEL METHOD
for iter = 1:maxIt
    %%% CREATE VELOCITY FROM STREAM FUNCTION
    u(i, j) = (St(i, jp) - St(i, jm)) / (2*h);
    v(i, j) = (-St(ip, j) + St(im, j)) / (2*h);

     T(1, 1:Ny) = T(2, 1:Ny);
    T(Nx, 1:Ny) = T(Nx-1, 1:Ny);
    T(1:Nx, 1) = 300; 
    T(1:Nx, Ny) = 330; 
    %%% TEMPERATURE EQUATION (ENERGY EQUATION)
    % Temporal discretization
    Tn = T;
    T(i,j) =(k * (Tn(ip, j) + Tn(im, j) + Tn(i, jp) + Tn(i, jm)) + ...
        max(u(i, j), 0) .* h .* Tn(im, j) + max(-u(i, j), 0) .* h .* Tn(ip, j) + ...
        max(v(i, j), 0) .* h .* Tn(i, jm) + max(-v(i, j), 0) .* h .* Tn(i, jp)) ./ ...
        (max(u(i, j), 0) .* h + max(-u(i, j), 0) .* h + max(v(i, j), 0) .* h + max(-v(i, j), 0) .* h + 4 * k);

 

    %%% SAVING THE VALUE OF VORTICITY AT THE LAST STEP
    Vop = Vo;

    %%% CREATE BOUNDARY CONDITIONS
    Vo(1:Nx, Ny) = -2 * St(1:Nx, Ny-1) / (h^2) - Wall_Velocity * 2 / h; % Top
    Vo(1:Nx, 1) = -2 * St(1:Nx, 2) / (h^2); % Bottom
    Vo(1, 1:Ny) = -2 * St(2, 1:Ny) / (h^2); % Left
    Vo(Nx, 1:Ny) = -2 * St(Nx-1, 1:Ny) / (h^2); % Right

    %%% PARTIALLY SOLVE VORTICITY TRANSPORT EQUATION
    Vo(i, j) = (mu / rho * (Vop(ip, j) + Vop(im, j) + Vop(i, jp) + Vop(i, jm)) + ...
        max(u(i, j), 0) .* h .* Vop(im, j) + max(-u(i, j), 0) .* h .* Vop(ip, j) + ...
        max(v(i, j), 0) .* h .* Vop(i, jm) + max(-v(i, j), 0) .* h .* Vop(i, jp)) ./ ...
        (max(u(i, j), 0) .* h + max(-u(i, j), 0) .* h + max(v(i, j), 0) .* h + max(-v(i, j), 0) .* h + 4 * mu / rho);

    %%% PARTIALLY SOLVE ELLIPTIC VORTICITY EQUATION FOR STREAM FUNCTION
    St(i, j) = (Vo(i, j) * h^2 + St(ip, j) + St(i, jp) + St(i, jm) + St(im, j)) / 4;


    %%% CHECK FOR CONVERGENCE
    if iter > 10
        error = max(max(Vo - Vop));
        if error < maxe
            break;
        end
    end
end

% Plot temperature contour
contourf(x, y, T', 20, 'LineWidth', 1.5);
colorbar;
colormap('jet');
xlabel('X');
ylabel('Y');
title('Temperature Contour');

