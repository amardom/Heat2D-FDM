% % % % % % % % % % % % % % % % % % %
%     HEAT EQUATION 1D UNSTEADY     %
%                                   %
%  dT/dt = alpha*dT/dx              %
%                                   %
%  T(x,0) = 0                       %
%  T(0,t) = 200                     %
%  dT/dx(L,t) = 300                 %
%                                   %
%  Finite Difference Method         %
%  - Time: central differences.     %
%  - Space: central differences.    %
%                                   %
%  A. Martínez                      %
% % % % % % % % % % % % % % % % % % %

clear;

% Initialize variables.
L = 1;
Nx = 100;
alpha = 1;
n_timestamps = 80000;
dx = L/Nx;
dt = 0.00001;
T = zeros(n_timestamps,Nx);
T_old = zeros(1,Nx);
T_old(1) = 200;
T_old(Nx) = 300;

% Stability condition.
VonNeu = alpha*dt/dx^2;
fprintf('\n ## VonNeu: %2.4f \n\n', VonNeu);
pause(2);
if (VonNeu > 0.5)
    fprintf(' ## Stopping program: VonNeu too high.\n');
    return;
end

for k = 1:n_timestamps
    
    for j = 2:Nx-1
        T(k,j) = T_old(j) + alpha*(T_old(j+1)/dx^2 - (2*T_old(j))/dx^2 + T_old(j-1)/dx^2)*dt + 300*dt;
    end
    
    T(k,1) = 300;
    T(k,Nx) = 400;
    
    T_old = T(k,:);
    
    disp(k);
end

% Plot.
x = linspace(0,L,Nx);
for k = 1:50:n_timestamps
    time = k*dt;
    plot(x,T(k,:));
    title(sprintf('%.4f s (timestamp: %i)',time,k));
    axis([0 L 280 420])
    pause(0.00001);
end