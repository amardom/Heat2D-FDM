% % % % % % % % % % % % % % % % % % %
%     WAVE EQUATION 1D UNSTEADY     %
%                                   %
%  d2u/dt2 = c^2*(d2u/dx2)          %
%                                   %
%  u(x,0) = 0                       %
%  u(0,t) = 200                     %
%  dT/dx(L,t) = 300                 %
%                                   %
%  Finite Difference Method         %
%  - Time: central differences      %
%  - Space: central differences     %
%                                   %
%  A. Martínez                      %
% % % % % % % % % % % % % % % % % % %

clear;

% Initialize variables.
L = 1;
Nx = 500;
c = 50;
n_timestamps = 4000;
dx = L/Nx;
dt = 0.00001;
u_p = zeros(Nx); % "p" = past, "n" = now, and "f" = future. 
u_n = zeros(Nx);
u_f = zeros(n_timestamps,Nx);
n_cycles = 1;

% CFL condition.
CFL = c*dt/dx;
fprintf('\n ## CFL: %2.4f \n\n', CFL);
if (CFL > 0.5)
    fprintf(' ## Stopping program: CFL too high.\n');
    return;
end
pause(2);

% Boundary conditions (Dirichlet)
u_f(:,1) = 0;
u_f(:,Nx) = 0;

% Compute u_d.
for k = 1:n_timestamps
    
    for j = 2:Nx-1
        u_f(k,j) = -u_p(j)...
                  + 2*u_n(j)...
                  + c^2*(u_n(j+1)/dx^2 - (2*u_n(j))/dx^2 + u_n(j-1)/dx^2)*dt^2;
    end
    
    if (0.1*k <= 2*pi*n_cycles)
        u_f(k,1) = sin(0.1*k);
    else
        u_f(k,1) = 0;
    end
    
    u_p = u_n;
    u_n = u_f(k,:);
    
    disp(k);
end

% Plot.
x = linspace(0,L,Nx);
for k = 1:10:n_timestamps
    time = k*dt;
    plot(x,u_f(k,:));
    title(sprintf('%.4f s (timestamp: %i)',time,k));
    axis([0 L -2 2])
    pause(0.05);
end
