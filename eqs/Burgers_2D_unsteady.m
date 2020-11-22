% % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
%               BURGERS EQUATION 2D UNSTEADY              % 
%                                                         %
%  du/dt + u*(du/dx) + u*(du/dy) = nu*(d2u/dx2 + d2u/dy2) %
%                                                         %
%  u(x,0) = 0                                             %
%  u(0,t) = 10                                            %
%  dT/dx(L,t) = 10                                        %
%                                                         %
%  Finite Difference Method                               %
%  - Time: forward differences.                           %
%  - Space: central differences.                          %
%                                                         %
%  A. Martínez                                            %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

clear;

% Initialize variables.
Lx = 1;
Ly = 1;
Nx = 100;
Ny = 100;
nu = 1;
n_timestamps = 1500;
dx = Lx/Nx;
dy = Ly/Ny;
dt = 0.00001;
u = zeros(n_timestamps,Ny,Nx);
u_old = zeros(Ny,Nx);
u_initial = 100;

% CFL condition.
CFL = u_initial*dt/dx;
fprintf('\n ## CFL: %2.4f \n\n', CFL);
if (CFL > 0.1)
    fprintf(' ## Stopping program: CFL too high.\n');
    return;
end
pause(2);

% Initial condition.
for j = (Nx/2-Nx*0.05):(Nx/2+Nx*0.05)
    for i = (Ny/2-Ny*0.05):(Ny/2+Ny*0.05)
        u_old(i,j) = u_initial;
    end
end

% BCs.
u(:,1,:) = 0;
u(:,Ny,:) = 0;
u(:,:,1) = 0;
u(:,:,Nx) = 0;

% Compute u.
for k = 1:n_timestamps
    
    for j = 2:Nx-1
        for i = 2:Ny-1
            
            u(k,i,j) = u_old(i,j)...
                   - u_old(i,j)*(u_old(i,j+1)/(2*dx) - u_old(i,j-1)/(2*dx))*dt...
                   - u_old(i,j)*(u_old(i+1,j)/(2*dy) - u_old(i-1,j)/(2*dy))*dt...
                   + nu*(u_old(i,j+1)/dx^2 - (2*u_old(i,j))/dx^2 + u_old(i,j-1)/dx^2)*dt...
                   + nu*(u_old(i+1,j)/dy^2 - (2*u_old(i,j))/dy^2 + u_old(i-1,j)/dy^2)*dt;
           
        end
    end
    
    u_old(:,:) = u(k,:,:);
    
    disp(k);
    
end

% Plot.
x = linspace(0,Lx,Nx);
y = linspace(0,Ly,Ny);
[X,Y] = meshgrid(y,x);

for k = 1:1:n_timestamps
    time = k*dt;
    
    single_snapshot = u(k,:,:);
    single_snapshot = squeeze(single_snapshot);
    single_snapshot = flipud(single_snapshot);
    
    figure(1);
    s = surf(X,Y,single_snapshot);
    title(sprintf('%.4f s (timestamp: %i)',time,k));
    colorbar; colormap jet;
    s.EdgeColor = 'interp';
    axis([0 1 0 1 -1 100]);
    caxis([-0.1 0.1]);
    pause(0.001);
end