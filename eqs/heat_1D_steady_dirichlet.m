% % % % % % % % % % % % % % % % % %
%     HEAT EQUATION 1D STEADY     %
%                                 %
%  0 = alpha*dT/dx + S(x)         %
%                                 %
%  T(x,0) = 0                     %
%  T(0,t) = 300                   %
%  T(L,t) = 400                   %
%                                 %
%  Finite Difference Method       %
%  - Space: central differences   %
%                                 %
%  A. Martínez                    %
% % % % % % % % % % % % % % % % % %

clear;

% Initialize variables.
L = 1;
Nx = 100;
alpha = 1;
dx = L/Nx;
rhs = zeros(Nx,1);
A = zeros(Nx);

% We construct the laplacian operator.
node = 1;
for j = 1:Nx
    
    if (j == 1 || j == Nx)
        node = node+1;
        continue;
    end
    
    A(node,node) = alpha*2/dx^2;
    A(node,node+1) = alpha*(-1)/dx^2;
    A(node,node-1) = alpha*(-1)/dx^2;
    
    node = node+1;
    
end

% Source term.
for j = 2:Nx-1
    rhs(j) = 300;
end

% Dirichlet BC.
A(1,1) = 1;
rhs(1) = 300;
A(Nx,Nx) = 1;
rhs(Nx) = 400;

% Solve.
T = A\rhs;

% Plot.
x = linspace(0,L,Nx);
figure(2);
plot(x,T');