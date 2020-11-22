clear;

% Input -START-
Lx = 0.1;
Ly = 0.1;
Nx = 4; 
Ny = 4; 

n_timestamps = 1000;
dt = 0.000005;  

u_TOP = 6.5; 
v_TOP = 0;

u_BOTTOM = 0;
v_BOTTOM = 0;

u_LEFT = 0;
v_LEFT = 0;

u_RIGHT = 0;
v_RIGHT = 0;

mu = 0.01;
rho = 1;
% Input -END-

dx = Lx/Nx;
dy = Ly/Ny;
dxi = 1/dx;
dyi = 1/dy;

nu = mu/rho;

reynolds = u_TOP*Lx/nu;
courant = dt*u_TOP/dx;

fprintf('\n\n ## Reynolds: %2.1f', reynolds);
fprintf('\n ## Courant: %2.4f \n\n', courant);
pause(1);

if (courant > 0.1)
    return;
end

u = zeros(Ny-1, Nx-1);
v = zeros(Ny-1, Nx-1);
u_star = zeros(Ny-1, Nx-1);
v_star = zeros(Ny-1, Nx-1);

rhs = zeros((Ny-1)*(Nx-1), 1);
A = zeros((Ny-1)*(Nx-1), (Ny-1)*(Nx-1));
p = zeros(Ny-1, Nx-1);

u_2plot = zeros(n_timestamps, Ny-1, Nx-1);
v_2plot = zeros(n_timestamps, Ny-1, Nx-1);
p_2plot = zeros(n_timestamps, Ny-1, Nx-1);
      
% Initialize coefficient matrix of the pressure equation.
node = 1;
for i = 1:Ny-1
    for j = 1:Nx-1

        if (i == 1)
            if(j == 1 || j == Nx-1) %Corner.
                A(node, node) = -2/dx^2 -2/dy^2;
                node = node+1;
                continue;
            end
            A(node, node) = -2/dx^2 -2/dy^2;
            A(node, node + 1) = 1/dx^2;
            A(node, node - 1) = 1/dx^2;
            A(node, node + (Nx-1)) = 2/dy^2;
            node = node+1;
            continue;
        end
        
        if (i == Ny-1)
            if(j == 1 || j == Nx-1) %Corners.
                A(node, node) = -2/dx^2 -2/dy^2;
                node = node+1;
                continue;
            end
            A(node, node) = -2/dx^2 -2/dy^2;
            A(node, node + 1) = 1/dx^2;
            A(node, node - 1) = 1/dx^2;
            A(node, node - (Nx-1)) = 2/dy^2;
            node = node+1;
            continue;
        end
        
        if (j == 1)
            if(i == 1 || i == Ny-1) %Corners.
                A(node, node) = -2/dx^2 -2/dy^2;
                node = node+1;
                continue;
            end
            A(node, node) = -2/dx^2 -2/dy^2;
            A(node, node + (Nx-1)) = 1/dy^2;
            A(node, node - (Nx-1)) = 1/dy^2;
            A(node, node + 1) = 2/dx^2;
            node = node+1;
            continue;
        end
        
        if (j == Nx-1)
            if(i == 1 || i == Ny-1) %Corners.
                A(node, node) = -2/dx^2 -2/dy^2;
                node = node+1;
                continue;
            end
            A(node, node) = -2/dx^2 -2/dy^2;
            A(node, node + (Nx-1)) = 1/dy^2;
            A(node, node - (Nx-1)) = 1/dy^2;
            A(node, node - 1) = 2/dx^2;
            node = node+1;
            continue;
        end
        
        A(node, node) = -2/dx^2 -2/dy^2;
        A(node, node + 1) = 1/dx^2;
        A(node, node - 1) = 1/dx^2;
        A(node, node + (Nx-1)) = 1/dy^2;
        A(node, node - (Nx-1)) = 1/dy^2;
        
        node = node+1;
        
    end
end

deter1 = det(A);
A(1,:) = 0;
A(1,1) = 1;
deter2 = det(A);

for k = 1:n_timestamps
    
    % BCs for the velocity.
    itop = 1;
    ibottom = Ny-1;
    jleft = 1;
    jright = Nx-1;
    
    for j = 1:Nx-1
        
        %u(itop,j) = u(itop+1,j) - 2*(u(itop+1,j) - u_TOP);
        u(itop,j) = -u(itop+1,j) + 2*u_TOP;
        v(itop,j) = v_TOP;
        
        %u(ibottom,j) = u(ibottom-1,j) - 2*(u(ibottom-1,j) - u_BOTTOM);
        u(ibottom,j) = -u(ibottom-1,j) + 2*u_BOTTOM;
        v(ibottom,j) = v_BOTTOM;
        
    end
    
    for i = 1:Ny-1
        
        u(i,jleft) = u_LEFT;
        %v(i,jleft) = v(i,jleft+1) - 2*(v(i,jleft+1) - v_LEFT);
        v(i,jleft) = -v(i,jleft+1) + 2*v_LEFT;
        
        u(i,jright) = u_RIGHT;
        %v(i,jright) = v(i,jright-1) - 2*(v(i,jright-1) - v_RIGHT);
        v(i,jright) = -v(i,jright-1) + 2*v_RIGHT;
        
    end
    
    % Compute u_star, and v_star;
    for j = 2:Nx-2
        for i = 2:Ny-2

            v_this_point = 0.25*(v(i-1,j) + v(i-1,j+1) + v(i,j) + v(i,j+1));
            u_star(i,j) = u(i,j) + dt*...
                        ( nu*(u(i-1,j) - 2*u(i,j) + u(i+1,j))*dxi^2 ...
                        + nu*(u(i,j-1) - 2*u(i,j) + u(i,j+1))*dyi^2 ...
                        - u(i,j)*(u(i+1,j) - u(i-1,j))*0.5*dxi ...
                        - v_this_point*(u(i,j+1) - u(i,j-1))*0.5*dyi);
                    
            u_this_point = 0.25*(u(i,j-1) + u(i,j) + u(i+1,j-1) + u(i+1,j));
            v_star(i,j) = v(i,j) + dt* ...
                        ( nu*(v(i-1,j) - 2*v(i,j) + v(i+1,j))*dxi^2 ...
                        + nu*(v(i,j-1) - 2*v(i,j) + v(i,j+1))*dyi^2 ...
                        - u_this_point*(v(i+1,j) - v(i-1,j))*0.5*dxi ...
                        - v(i,j)*(v(i,j+1) - v(i,j-1))*0.5*dyi);

        end
    end

    % rhs pressure equation.
    node = 1;
    for j = 1:Nx-1
        for i = 1:Ny-1
           
            if (i == 1 || i == Ny-1 || j == 1 || j == Nx-1)
                node = node+1;
                continue;
            end
            
            rhs(node) = -rho/dt * ...
                        ((u_star(i+1,j) - u_star(i,j)) * dxi ...
                        +(v_star(i,j+1) - v_star(i,j)) * dyi);
            
            node = node+1;
            
        end
    end
    
    % Solve pressure equation.
    p_vector = A\rhs;

    % Pressure vector to matrix.
    node = 1;
    for j = 1:Nx-1
        for i = 1:Ny-1
            
            p(i,j) = p_vector(node);
            node = node + 1;
            
        end
    end

    % Correct u, and v.
    for j = 2:Nx-2
        for i = 2:Ny-2

            u(i,j) = u_star(i,j) - dt/rho * (p(i,j) - p(i-1,j)) * dxi;
            v(i,j) = v_star(i,j) - dt/rho * (p(i,j) - p(i,j-1)) * dyi;

        end
    end
    
    % Store variables.
    u_2plot(k,:,:) = u;
    v_2plot(k,:,:) = v;
    p_2plot(k,:,:) = p;
    disp(k);
    
end

% Plot.
x = linspace(0,Lx,Nx-1);
y = linspace(0,Ly,Ny-1);
[X,Y] = meshgrid(x,y);

for k = 1:1:n_timestamps
    time = k*dt;
    v_mag = sqrt(u_2plot(k,:,:).^2 + v_2plot(k,:,:).^2);
    
    single_snapshot = u_2plot(k,:,:);
    single_snapshot = squeeze(single_snapshot);
    single_snapshot = flipud(single_snapshot);
    
    figure(2);
    contourf(X,Y,single_snapshot,10);
    title(sprintf('%.4f s (timestamp: %i)',time,k));
    colorbar; colormap jet;
    set(gca,'DataAspectRatio',[1 1 1]);
    set(gcf,'position',[1000,100,850,850]);
    pause(0.01);
end