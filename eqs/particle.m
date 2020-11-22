clear

n_timestamps = 8000;
g = 9.8;
x = zeros(n_timestamps,1);
y = zeros(n_timestamps,1);
dt = 0.001;
x0 = 0;
y0 = 0.5;
vx0 = 15;
vy0 = -5;
t = 0;

BOTTOM = 0;
TOP = 1;
RIGHT = 1;
LEFT = -1;

for k = 1:n_timestamps
    
    if (k < 500)
        x(k) = x0; 
        y(k) = y0;
    else
        t = t + dt;
        x(k) = x0 + vx0*t;
        y(k) = y0 + vy0*t - g*(t/2)^2;
    end
    
    if (y(k) < BOTTOM)
        
        vx0 = -((x(k-1) - x(k))/dt);
        vy0 = ((y(k-1) - y(k))/dt);
        x0 = x(k);
        y0 = y(k);
        t = 0;
        
    end
    
    if (y(k) > TOP)
        
        vx0 = -((x(k-1) - x(k))/dt);
        vy0 = ((y(k-1) - y(k))/dt);
        x0 = x(k);
        y0 = y(k);
        t = 0;
        
    end
    
    if (x(k) > RIGHT)

        vx0 = ((x(k-1) - x(k))/dt);
        vy0 = -((y(k-1) - y(k))/dt);
        x0 = x(k);
        y0 = y(k);
        t = 0;

    end

    if (x(k) < LEFT)

        vx0 = ((x(k-1) - x(k))/dt);
        vy0 = -((y(k-1) - y(k))/dt);
        x0 = x(k);
        y0 = y(k);
        t = 0;

    end
    
end


for k = 1:2:n_timestamps
    plot(x(k),y(k),'.','MarkerSize',30);
    axis([-1 1 0 1]);
    pause(0.001);
end