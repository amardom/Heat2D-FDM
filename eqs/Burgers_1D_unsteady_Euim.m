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
%  - Time: backward differences.      %
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
u_iter = zeros(1,Nx);
u_initial = 100;
A = zeros(Nx,Nx);
rhs = zeros(Nx,1);
semi_throat = Nx*0.4;
n_iter = 5;

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

% BCs.
A(1,1) = 1;
A(Nx,Nx) = 1;

% Compute u.
for k = 1:n_timestamps
    
    for iter = 1:n_iter % Picard.

        for j = 2:Nx-1

            A(j,j) = 1/dt + 2*nu*(1/dx^2); 
            A(j,j+1) = u_iter(j)/(2*dx) - nu*(1/dx^2);
            A(j,j-1) = - u_iter(j)/(2*dx) - nu*(1/dx^2);

        end

        rhs = u_old'/dt;

        u(k,:) = A\rhs;
        u_iter = u(k,:);
        
        fprintf(' ##Timestamp:%d  norm:%8.12f \n', k, norm(u_iter));
        
    end
    
    disp(' ');
    
    u_old = u_iter;
    
    u_old(1) = 0;
    u_old(Nx) = 0;
    
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
