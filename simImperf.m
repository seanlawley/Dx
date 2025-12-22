function output = simImperf(a,alpha,b0,b1,kappa)
% parameters
Lx = 1; Ly = 1; Lz = 1;
dt_min = 1e-8;
dt_max = 5e-3;
beta   = 0.1;

% D(x)=b0+(b1-b0)*x

diskRadius2 = a^2;
yc = Ly/2;
zc = Lz/2;

% Initial position
x = Lx/2;
y = Ly/2;
z = Lz/2;

t = 0;

% simulation loop
while true

    % Distance to absorbing region (approximate):
    % radial distance in y-z
    dy = y - yc; dz = z - zc;
    rho = sqrt(dy*dy + dz*dz);
    % distance from current point to disk boundary (0 when inside)
    dist_to_disk = sqrt( min(x, Lx - x)^2 + max(0, rho - a)^2 );
    if dist_to_disk <= 0
        % we're (practically) at the disk boundary: use smallest dt
        dt = dt_min;
    else
        % choose dt so RMS step ~ beta * dist_to_disk: dt = (beta*dist)^2/(2D)
        dt_choice = (beta * dist_to_disk)^2 / (2*(b0+(b1-b0)*x));
        % clamp
        dt = min(max(dt_choice, dt_min), dt_max);
    end

    % Vectorized 3D step
    r = sqrt((b0+(b1-b0)*x)*2*dt) * randn(1,3);
    x = x + alpha * (b1-b0) * dt + r(1);
    y = y + r(2);
    z = z + r(3);

    % Reflect y and z
    if y < 0 || y > Ly
        y = abs(mod(y, 2*Ly));
        if y > Ly, y = 2*Ly - y; end
    end

    if z < 0 || z > Lz
        z = abs(mod(z, 2*Lz));
        if z > Lz, z = 2*Lz - z; end
    end

    % absorbing disks
    dy = y - yc;
    dz = z - zc;
    r2 = dy*dy + dz*dz;

    if r2 <= diskRadius2   
        % hit left?
        if x <= 1e-12 && binornd(1,kappa*sqrt(dt*pi/b0))==1
            T = t;
            exitRight = 0;
            output = [T exitRight];
            return
        end
        % hit right?
        if Lx - x <= 1e-12 && binornd(1,kappa*sqrt(dt*pi/b1))==1
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
