% % % % % % % % % % % % % % % % % % % % % %
%     WAVE EQUATION 2D UNSTEADY           %
%                                         %
%  d2u/dt2 = c^2*(d2u/dx2 + d2u/dy2)      %
%                                         %
%  u(x,y,0) = 0                           %
%  u(0,y,t) = u(Lx,y,t) = 0               %
%  u(x,0,t) = u(x,Ly,t) = 0               %
%  u(Lx/2,Ly/2,t) = 2*sin(t) for 0<t<Tst  %
%                                         %
%  Finite Difference Method               %
%  - Time: forward differences            %
%  - Space: central differences           %
%                                         %
%  A. Martínez                            %
% % % % % % % % % % % % % % % % % % % % % %

clear;

% Initialize variables.
Lx = 1;
Ly = 1;
Nx = 500;
Ny = 500;
c = 50;
n_timestamps = 2000;
dx = Lx/Nx;
dy = Ly/Ny;
dt = 0.00001;
u_p = zeros(Ny,Nx); % "p" = past, "n" = now, and "f" = future. 
u_n = zeros(Ny,Nx);
u_f = zeros(n_timestamps,Ny,Nx);
n_cycles = 5;

% CFL condition.
CFL = c*dt/dx;
fprintf('\n ## CFL: %2.4f \n\n', CFL);
if (CFL > 0.5)
    fprintf(' ## Stopping program: CFL too high.\n');
    return;
end
pause(2);

% Boundary conditions (Dirichlet).
u_f(:,1,:) = 0;
u_f(:,Ny,:) = 0;
u_f(:,:,1) = 0;
u_f(:,1,Nx) = 0;

% Compute u_f.
for k = 1:n_timestamps
    
    for j = 2:Nx-1
        for i = 2:Ny-1
            
            u_f(k,i,j) = -u_p(i,j)...
                             + 2*u_n(i,j)...
                             + c^2*(u_n(i,j+1)/dx^2 - (2*u_n(i,j))/dx^2 + u_n(i,j-1)/dx^2)*dt^2 ...
                             + c^2*(u_n(i+1,j)/dy^2 - (2*u_n(i,j))/dy^2 + u_n(i-1,j)/dy^2)*dt^2;
                     
        end
    end
    
    if (0.1*k <= 2*pi*n_cycles)
        u_f(k,Nx/2,Ny/2) = 2*sin(0.1*k);
    end
    
    u_p = u_n;
    u_n = squeeze(u_f(k,:,:));
    
    disp(k);
    
end

% Plot.
x = linspace(0,Lx,Nx);
y = linspace(0,Ly,Ny);
[X,Y] = meshgrid(y,x);

for k = 1:20:n_timestamps
    time = k*dt;
    
    single_snapshot = u_f(k,:,:);
    single_snapshot = squeeze(single_snapshot);
    single_snapshot = flipud(single_snapshot);
    
    figure(1);
    s = surf(X,Y,single_snapshot);
    title(sprintf('%.4f s (timestamp: %i)',time,k));
    colorbar; colormap jet;
    s.EdgeColor = 'interp';
    axis([0 1 0 1 -1 1]);
    caxis([-0.1 0.1]);
    pause(0.001);
end
