
%% VORTICITY/STREAM FUNCTION LID DRIVEN CAVITY FLOW SOLVER JOE MOLVAR
clear; 
clc;
close all
%%% GIVENS
Nx = 52; %Nodes 
L = 1; %Domain Size
rho = 1;% Density
mu = 0.01; %Dynamic Viscosity
Re=100; %Reynolds
Wall_Velocity =1; %Velocity
dt = 0.001;%Time Step
maxIt = 50000; %Max iter
maxe = 1e-7; %Max error
%%% SETUP 1D GRID
Ny = Nx; h=L/(Nx-1); x = 0:h:L; y = 0:h:L;
im = 1:Nx-2; i = 2:Nx-1; ip = 3:Nx; jm = 1:Ny-2; j = 2:Ny-1; jp = 3:Ny; I=1:Nx; J=1:Ny;
%%% PRELOCATE MATRIXES
Vo = zeros(Nx,Ny); St = Vo; Stp = Vo; Vop = Vo; u = Vo; v = Vo; P = Vo; epsilon = 0;
%%% VELOCITY ON THE UPPER WALL(NO SLIP CONDITION)
u(2:Nx-1,Ny) = Wall_Velocity;
%%% SOLVE LOOP SIMILAR TO GAUSS-SIEDEL METHOD
for iter = 1:maxIt
%%% CREATE BOUNDARY CONDITIONS
Vo(1:Nx,Ny) = -2*St(1:Nx,Ny-1)/(h^2) - Wall_Velocity*2/h; % Top
Vo(1:Nx,1) = -2*St(1:Nx,2) /(h^2); % Bottom
Vo(1,1:Ny) = -2*St(2,1:Ny) /(h^2); % Left
Vo(Nx,1:Ny) = -2*St(Nx-1,1:Ny)/(h^2); % Right
%%% PARTIALLY SOLVE VORTICITY TRANSPORT EQUATION
Vop = Vo;
Stp = St;
Vo(i,j) = Vop(i,j) + ...
(-1*(St(i,jp)-St(i,jm))/(2*h) .* (Vop(ip,i)-Vop(im,j))/(2*h)+...
(St(ip,j)-St(im,j))/(2*h) .* (Vop(i,jp)-Vop(i,jm))/(2*h)+...
mu/rho*(Vop(ip,j)+Vop(im,j)-4*Vop(i,j)+Vop(i,jp)+Vop(i,jm))/(h^2))*dt;
%%% PARTIALLY SOLVE ELLIPTICAL VORTICITY EQUATION FOR STREAM FUNCTION
St(i,j) = (Vo(i,j)*h^2 + St(ip,j) + St(i,jp) + St(i,jm) + St(im,j))/4;
%%% CREATE VELOCITY FROM STREAM FUNCTION
u(i,j) = (St(i,jp)-St(i,jm))/(2*h); v(i,j) = (-St(ip,j)+St(im,j))/(2*h);
%%% PRESSURE TOP BOUNDARY CONDITION
P(I,Ny) = 1/3*(4*P(I,Ny-1)-P(I,Ny-2))+(2*mu)/(3*h)*(-5*v(I,Ny-1)+4*v(I,Ny-2)-v(I,Ny-3));
%%% PRESSURE BOTTOM BOUNDARY CONDITION
P(I,1) = 1/3*(4*P(I,2)-P(I,3))-(2*mu)/(3*h)*(-5*v(I,2)+4*v(I,3)-v(I,4));
%%% PRESSURE RIGHT BOUNDARY CONDITION
P(Nx,J) = 1/3*(4*P(Nx-1,J)-P(Nx-2,J))+(2*mu)/(3*h)*(-5*u(Nx-1,J)+4*u(Nx-2,J)-u(Nx-3,J));
%%% PRESSURE LEFT BOUNDARY CONDITION
P(1,J) = 1/3*(4*P(2,J)-P(3,J))-(2*mu)/(3*h)*(-5*u(2,J)+4*u(3,J)-u(4,J));
%%% PRESSURE FOR INNER NODES
P(i,j) = 0.25*(P(ip,j)+P(im,j)+P(i,jp)+P(i,jm))-rho/2*(1/(h^2)*(St(ip,j)-2*St(i,j)+....
St(im,j)).*(St(i,jp)-2*St(i,j)+St(i,jm))-1/(16*h^2)*(St(ip,jp)-St(ip,jm)-St(im,jp)+St(im,jm)).^2);
%%% CALCULATING EPSILON
epsilon(1,iter) = max(abs(St-Stp),[],'all');
%%% CHECK FOR CONVERGENCE
if iter > 10
error = max(max(Vo - Vop));
if error < maxe
break;
end
end
end
%%% PLOTS
cm = hsv(ceil(100/0.7)); cm = flipud(cm(1:100,:));
figure(1); contourf(x,y,u',23,'LineColor','none');
title('U-velocity'); xlabel('x-location'); ylabel('y-location')
axis('equal',[0 L 0 L]); colormap(cm); colorbar('westoutside');
figure(2); plot(y,u(round(Ny/2),:));
title('Centerline x-direction velocity');
xlabel('y/L'); ylabel('u/U'); axis('square'); xlim([0 L]); grid on
N = 1000; xstart = max(x)*rand(N,1); ystart = max(y)*rand(N,1);
[X,Y] = meshgrid(x,y);
figure(3); h=streamline(X,Y,u',v',xstart,ystart,[0.1, 200]);
title('Stream Function'); xlabel('x-location'); ylabel('y-location')
axis('equal',[0 L 0 L]); set(h,'color','k')
figure(4); 
contourf(x,y,P',20,'LineColor','none');
title('Pressure'); xlabel('x-location'); ylabel('y-location'); axis('equal',[0 L 0 L]); colormap(cm); colorbar('westoutside');
figure(5); 
subplot(2,2,1);
plot(x,Vo(:,1),'LineWidth',2);
title('Vorticity on the bottom wall'); xlabel('x-location'); ylabel('Vorticity'); ylim([-80 10]);
grid on
subplot(2,2,2);
plot(x,Vo(:,Ny),'LineWidth',2);
title('Vorticity on the top wall'); xlabel('x-location'); ylabel('Vorticity'); ylim([-80 10]);
grid on
subplot(2,2,3);
plot(y,Vo(1,:),'LineWidth',2);
title('Vorticity on the left wall'); xlabel('y-location'); ylabel('Vorticity'); ylim([-2 40]);
grid on
subplot(2,2,4);
plot(y,Vo(Nx,:),'LineWidth',2);
title('Vorticity on the right wall'); xlabel('y-location'); ylabel('Vorticity'); ylim([-2 40]);
grid on
figure(6); 
contourf(x,y,Vop',20,'LineColor','none');
title('Vorticity'); xlabel('x-location'); ylabel('y-location'); axis('equal',[0 L 0 L]); colormap(cm); colorbar('westoutside');


