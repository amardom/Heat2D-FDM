clear

%----
n_timestamps = 30000;
dt = 0.00001;
g = 9.8;
n_balls = 201;
radius = 0.005;
m = ones(n_balls,1);
x = zeros(n_timestamps,n_balls,1);
y = zeros(n_timestamps,n_balls,1);
vx = zeros(n_timestamps,n_balls,1);
vy = zeros(n_timestamps,n_balls,1);
x0 = zeros(n_balls,1);
y0 = zeros(n_balls,1);
vx0 = zeros(n_balls,1);
vy0 = zeros(n_balls,1);
t = zeros(n_balls,1);
BOTTOM = 0;
TOP = 1;
RIGHT = 1;
LEFT = -1;
x_ini = -1;
y_ini_top = 0.95;
y_ini_bottom = 0.05;
y_ini = y_ini_top;
%----

j = 0;
for i = 1:n_balls
    j = j + 1;
    x0(i) = x_ini + (2/11)*j;
    y0(i) = y_ini;
    if (mod(i,10) == 0)
        y_ini = y_ini - (y_ini_top-y_ini_bottom)/(fix((n_balls)/10));
        j = 0;
    end
    vx0(i) = ((mod(fix(rand(1)*10),4))+3)*((-1)^(mod(fix(rand(1)*10),2)));
    vy0(i) = ((mod(fix(rand(1)*10),4))+3)*((-1)^(mod(fix(rand(1)*10),2)));
end
m(n_balls) = 500;

%{
x0(1) = 0;
y0(1) = 0.1;
vx0(1) = 0;
vy0(1) = 0;

x0(2) = 0;
y0(2) = 0.1;
vx0(2) = 0;
vy0(2) = 0;
%}

for k = 1:n_timestamps
    disp(k);
    collision_mask = zeros(n_balls,1);
    
    for i = 1:n_balls

        if (k < 100)
            x(k,i) = x0(i); 
            y(k,i) = y0(i);
        else
            t(i) = t(i) + dt;

            x(k,i) = x0(i) + vx0(i)*t(i);
            y(k,i) = y0(i) + vy0(i)*t(i) - g*0.5*t(i)^2;
            vx(k,i) = vx0(i);
            vy(k,i) = vy0(i) - g*t(i);
        end
        
    end
    
    for i = 1:n_balls
        
        for ii = 1:n_balls
           
            if (i == ii || collision_mask(i) == 1 || collision_mask(ii) == 1)
                continue;
            end

            r1(1) = x(k,i);
            r1(2) = y(k,i);
            
            r2(1) = x(k,ii);
            r2(2) = y(k,ii);
            
            distance = (r1(1)-r2(1))^2 + (r1(2)-r2(2))^2;

            if (distance < (radius+radius)^2)
                x0(i) = x(k,i);
                y0(i) = y(k,i);
                x0(ii) = x(k,ii);
                y0(ii) = y(k,ii);
                
                vbe1(1) = vx(k,i);
                vbe1(2) = vy(k,i);
                
                vbe2(1) = vx(k,ii);
                vbe2(2) = vy(k,ii);
                
                if (abs(vbe1(1) - vbe2(1)) < 0.0001 && abs(vbe1(2) - vbe2(2)) < 0.0001)                
                    y0(ii) = y0(ii) + 0.01*radius;
                    disp('HERE');
                end
                
                vaf1 = vbe1 - (2*m(ii)/(m(i)+m(ii))) * (dot(vbe1-vbe2,r1-r2)/norm(r1-r2)^2) * (r1-r2);
                vaf2 = vbe2 - (2*m(i)/(m(i)+m(ii))) * (dot(vbe2-vbe1,r2-r1)/norm(r2-r1)^2) * (r2-r1);
                
                vx0(i) = vaf1(1);
                vy0(i) = vaf1(2);
                vx0(ii) = vaf2(1);
                vy0(ii) = vaf2(2);
                
                collision_mask(i) = 1;
                collision_mask(ii) = 1;
                t(i) = 0;
                t(ii) = 0;                
            end
            
        end
        
        if (y(k,i) < BOTTOM)       
            x0(i) = x(k,i);
            y0(i) = BOTTOM;
            vx0(i) = vx(k,i);
            vy0(i) = -vy(k,i);
            t(i) = 0;
        end

        if (y(k,i) > TOP)
            x0(i) = x(k,i);
            y0(i) = TOP;
            vx0(i) = vx(k,i);
            vy0(i) = -vy(k,i);
            t(i) = 0;
        end

        if (x(k,i) > RIGHT)
            x0(i) = RIGHT;
            y0(i) = y(k,i);
            vx0(i) = -vx(k,i);
            vy0(i) = vy(k,i);
            t(i) = 0;
        end

        if (x(k,i) < LEFT)
            x0(i) = LEFT;
            y0(i) = y(k,i);
            vx0(i) = -vx(k,i);
            vy0(i) = vy(k,i);
            t(i) = 0;
        end
        
    end
    
end

figure('units','normalized','outerposition',[0 0 1 1]);
for k = 1:100:n_timestamps
    time = k*dt;
    plot(x(k,1:n_balls-1),y(k,1:n_balls-1),'b.',x(k,n_balls),y(k,n_balls),'r.','MarkerSize',30);
    title(sprintf('%.4f s (timestamp: %i)',time,k));
    axis([-1 1 0 1]);
    pause(0.001);
end
