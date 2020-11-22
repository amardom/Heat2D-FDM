clear;

% Input -------------------- START --------------------
Nx = 30; %35 for Re = 400.
Ny = 30; %35 for Re = 400.
Lx = 0.1;
Ly = 0.1;

n_timestamps = 5000; %15000 for Re = 400.
dt = 0.00005;  %0.000007 for Re = 400.

ubot = 0;
vbot = 0;

utop = 6.5; %40 for Re = 400.
vtop = 0;

uleft = 0;
vleft = 0;

uright = 0;
vright = 0;

nu = 0.01;
rho = 1;
% Input -------------------- END --------------------

dx = Lx/Nx;
dy = Ly/Ny;
dxi = 1/dx;
dyi = 1/dy;

reynolds = utop*Lx/nu;
courant = dt*utop/dx;

fprintf('\n\n ######### Reynolds: %2.2f', reynolds);
fprintf('\n\n ######### Courant: %2.2f \n\n', courant);
pause(4);

if (courant > 0.1)
    return;
end

imin = 2; 
imax = imin + Nx - 1;
jmin = 2; 
jmax = jmin + Ny - 1;

mu = nu/rho;

u = zeros(imax+1, jmax+1);
v = zeros(imax+1, jmax+1);
us = zeros(imax+1, jmax+1);
vs = zeros(imax+1, jmax+1);

R = zeros(Nx*Ny,1);
L = zeros(Nx*Ny,Nx*Ny);

% Laplacian.
for j = 1:Ny
    for i = 1:Nx
        
        L(i+(j-1)*Nx, i+(j-1)*Nx) = 2*dxi^2 + 2*dyi^2;
        
        for ii = i-1:2:i+1
            if (ii > 0 && ii <= Nx) % Interior point
                L(i+(j-1)*Nx, ii+(j-1)*Nx) = -dxi^2;
            else % Neumann conditions on boundary
                L(i+(j-1)*Nx, i+(j-1)*Nx) = ...
                L(i+(j-1)*Nx, i+(j-1)*Nx) - dxi^2;
            end
        end
        
        for jj = j-1:2:j+1
            if (jj > 0 && jj <= Ny) % Interior point
                L(i+(j-1)*Nx, i+(jj-1)*Nx)= -dyi^2;
            else % Neumann conditions on boundary
                L(i+(j-1)*Nx, i+(j-1)*Nx) = ...
                L(i+(j-1)*Nx, i+(j-1)*Nx) - dyi^2;
            end
        end
        
    end
end

%L(1,:) = 0; % Pressure for first cell.
L(1,1) = 1;

for k = 1:n_timestamps
    
    %BC u.
    for i = imin:imax
        u(i,jmin-1) = u(i,jmin) - 2*(u(i,jmin) - ubot);
        u(i,jmax+1) = u(i,jmax) - 2*(u(i,jmax) - utop);
        v(i,jmin) = vbot;
        v(i,jmax) = vtop;
    end
    
    %BC v.
    for j = jmin:jmax
        v(imin-1,j) = v(imin,j) - 2*(v(imin,j) - vleft);
        v(imax+1,j) = v(imax,j) - 2*(v(imax,j) - vright);
        u(imin,j) = uleft;
        u(imax,j) = uright;
    end
    
    for j = jmin:jmax
        for i = imin+1:imax

            vhere = 0.25*(v(i-1,j) + v(i-1,j+1) + v(i,j) + v(i,j+1));

            us(i,j) = u(i,j) + dt* ...
            ( nu*(u(i-1,j) - 2*u(i,j) + u(i+1,j))*dxi^2 ...
            + nu*(u(i,j-1) - 2*u(i,j) + u(i,j+1))*dyi^2 ...
            - u(i,j)*(u(i+1,j) - u(i-1,j)) * 0.5*dxi ...
            - vhere*(u(i,j+1) - u(i,j-1)) * 0.5*dyi);

        end
    end

    for j = jmin+1:jmax
        for i = imin:imax

            uhere = 0.25*(u(i,j-1) + u(i,j) + u(i+1,j-1) + u(i+1,j));

            vs(i,j) = v(i,j) + dt* ...
            ( nu*(v(i-1,j) - 2*v(i,j) + v(i+1,j))*dxi^2 ...
            + nu*(v(i,j-1) - 2*v(i,j) + v(i,j+1))*dyi^2 ...
            - uhere*(v(i+1,j) - v(i-1,j)) * 0.5*dxi ...
            - v(i,j)*(v(i,j+1) - v(i,j-1)) * 0.5*dyi);
        end
    end

    %rhs P.
    n = 0;
    for j = jmin:jmax
        for i = imin:imax

            n = n + 1;
            R(n) = -rho/dt * ...
            ( (us(i+1,j) - us(i,j)) * dxi ...
             +(vs(i,j+1) - vs(i,j)) * dyi);

        end
    end

    %solve.
    pv=L\R;

    %vector representation to matrix.
    n = 0;
    p = zeros(imax,jmax);

    for j=jmin:jmax
        for i=imin:imax

            n = n + 1;
            p(i,j) = pv(n);

        end
    end

    %transform us,vs into u,v (correct it).
    for j = jmin:jmax
        for i = imin+1:imax

            u(i,j) = us(i,j) - dt/rho*(p(i,j) - p(i-1,j)) * dxi;

        end
    end

    for j = jmin+1:jmax
        for i = imin:imax

            v(i,j) = vs(i,j) - dt/rho*(p(i,j) - p(i,j-1)) * dyi;

        end
    end

    u_2plot(k,:,:) = u;
    v_2plot(k,:,:) = v;
    p_2plot(k,:,:) = p;
    disp(k);
end

% Plot.
x = linspace(0,Lx,imax+1);
y = linspace(0,Ly,jmax+1);
[X,Y] = meshgrid(y,x);

for k = 1:100:n_timestamps
    time = k*dt;
    v_mag = sqrt(u_2plot(k,:,:).^2 + v_2plot(k,:,:).^2);
    %Vorticity.
    [curlz, cav] = curl(X,Y,squeeze(u_2plot(k,:,:)),squeeze(v_2plot(k,:,:)));
    for i = 1:length(curlz(:,1))
        for j = 1:length(curlz(1,:))
            if (abs(curlz(i,j)) > 1000)
               curlz(i,j) = 0;
            end
        end
    end
    
    single_snapshot = v_mag;
    single_snapshot = squeeze(single_snapshot);
    single_snapshot = flipud(single_snapshot);
    
    figure(2);
    contourf(X,Y,single_snapshot,10);
    title(sprintf('%.4f s (timestamp: %i)',time,k));
    colorbar; colormap jet;
    set(gca,'DataAspectRatio',[1 1 1]);
    set(gcf,'position',[1000,100,850,850]);
    pause(0.0001);
end
