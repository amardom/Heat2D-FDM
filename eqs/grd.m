clear;

%-----
n_steps = 30000;
h = 0.05;
x = -2:h:2;
y = -3:h:3;
y = y';
Kgrd = 0.1;
max_diff = 1e-6;
%----

z = x.*exp(-x.^2-y.^2);
[grad_x, grad_y] = gradient(z);

a(1,1) = 0.65;
a(1,2) = -0.4;

for k = 2:n_steps
    
    disp(k);
    
    diff_old = 1e15;
    i_index = 0;
    for i = 1:length(x)
        diff = abs(x(i) - a(k-1,1));
        if (diff < diff_old)
            i_index = i;
            diff_old = diff;
        end
    end
    
    diff_old = 1e15;
    j_index = 0;
    for i = 1:length(y)
        diff = abs(y(i) - a(k-1,2));
        if (diff < diff_old)
            j_index = i;
            diff_old = diff;
        end
    end
    
    a(k,1) = a(k-1,1) - Kgrd*grad_x(j_index,i_index);
    a(k,2) = a(k-1,2) - Kgrd*grad_y(j_index,i_index);
    
    z_now = a(k,1)*exp(-a(k,1)^2 - a(k,2)^2);
    z_old = a(k-1,1)*exp(-a(k-1,1)^2 - a(k-1,2)^2);
    
    if (abs(z_now - z_old) < max_diff)
        break;
    end
    
end

figure('units','normalized','outerposition',[0 0 1 1]);
mesh(x,y,z);
hold on;

for k = 2:5:n_steps
    title(sprintf('GRADIENT DESCENT. step: %i',k));
    z = a(k,1)*exp(-a(k,1)^2 - a(k,2)^2);
    plot3(a(k,1),a(k,2),z,'r.','MarkerSize',20);
    pause(0.01);
end
