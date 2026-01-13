function output = simPerf2d(a,alpha,b0,b1)
Lx = 1; Ly = 1;
dt_min = min([.01*a^2/(2*max([b0,b1])),5e-3]);
dt_max = 5e-3;
beta   = 0.1;

diskRadius2 = a^2;
yc = Ly/2;

% Initial position
x = Lx/2;
y = Ly/2;

t = 0;

while true

    dy = y - yc;
    rho = abs(dy);
    dist_to_disk = sqrt( min(x, Lx - x)^2 + max(0, rho - a)^2 );
    if dist_to_disk <= 0
        dt = dt_min;
    else
        dt_choice = (beta * dist_to_disk)^2 / (2*(b0+(b1-b0)*x));
        dt = min(max(dt_choice, dt_min), dt_max);
    end

    r = sqrt((b0+(b1-b0)*x)*2*dt) * randn(1,2);
    x = x + alpha * (b1-b0) * dt + r(1);
    y = y + r(2);

    if y < 0 || y > Ly
        y = abs(mod(y, 2*Ly));
        if y > Ly, y = 2*Ly - y; end
    end

    dy = y - yc;
    r2 = dy*dy;

    if r2 <= diskRadius2  
        % hit left?
        if x <= 1e-12
            T = t;
            exitRight = 0;
            output = [T exitRight];
            return
        end
        % hit right?
        if Lx - x <= 1e-12
            T = t;
            exitRight = 1;
            output = [T exitRight];
            return
        end
    end

    % reflect x
    if x < 0 || x > Lx
        x = abs(mod(x, 2*Lx));
        if x > Lx, x = 2*Lx - x; end
    end

    t = t + dt;

end
end
