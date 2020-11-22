% % % % % % % % % % % % % % % % % % % %
%     BURGERS EQUATION 1D UNSTEADY    % 
%                                     %
%  du/dt + u*(du/dx) = nu*(d2u/dx2)   %
%                                     %
%  u(x,0) = 0                         %
%  u(0,t) = 10                        %
%  dT/dx(L,t) = 10                    %
%                                     %
%  Finite Difference Method           %
%  - Time: forward differences.       %
%  - Space: central differences.      %
%                                     %
%  A. Martínez                        %
% % % % % % % % % % % % % % % % % % % %

clear;

% Initialize variables.
L = 1;
Nx = 800;
nu = 0.1;
n_timestamps = 20000;
dx = L/Nx;
dt = 0.000001;
u = zeros(n_timestamps,Nx);
u_old = zeros(1,Nx);
u_initial = 100;
semi_throat = Nx*0.4;

% CFL condition.
CFL = u_initial*dt/dx;
fprintf('\n ## CFL: %2.4f \n\n', CFL);
if (CFL > 0.8)
    fprintf(' ## Stopping program: CFL too high.\n');
    return;
end
pause(2);

% Initial condition.
for j = (Nx/2-semi_throat):(Nx/2+semi_throat)
    u_old(j) = u_initial*exp(-(j-(Nx/2))^2/(2*semi_throat));
end

% Compute u.
for k = 1:n_timestamps
    
    for j = 2:Nx-1
        u(k,j) = u_old(j)...
               - u_old(j)*(u_old(j+1)/(2*dx) - u_old(j-1)/(2*dx))*dt...
               + nu*(u_old(j+1)/dx^2 - (2*u_old(j))/dx^2 + u_old(j-1)/dx^2)*dt;
    end
    
    u(k,1) = 0;
    u(k,Nx) = 0;
    
    u_old = u(k,:);
    
    disp(k);
    
end

% Plot.
x = linspace(0,L,Nx);
for k = 1:20:n_timestamps
    time = k*dt;
    plot(x,u(k,:));
    title(sprintf('%.4f s (timestamp: %i)',time,k));
    axis([0 L 0 u_initial*1.1])
    pause(0.001);
end
