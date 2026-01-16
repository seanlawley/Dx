function output = simTubes(params,alpha)

D1  = params.D1;
D2  = params.D2;
dt  = params.dtmin;
L0  = params.L0;
R0  = params.R0;
L1  = params.L1;
L2  = params.L2;
R1  = params.R1;
Tmax = params.Tmax;

% --------- Initial position in main cylinder ----------
x = 0; % initial x position along the main cylinder
r = R0*sqrt(rand); % initial radius
y = r;
z = 0;

t = 0; % initialize time to be zero

while t < Tmax

    % ---------- Determine region ----------
    % ----- and set diffusion coefficient ---
    if x >= -L0/2 && x <= L0/2
        region = "main";
        D = D1+(x+L0/2)*(D2-D1)/L0; % bulk diffusion
    elseif x < -L0/2 && r_new < R1
        region = "left";
        D = D1;
    elseif x > L0/2 && r_new < R1
        region = "right";
        D = D2;
    else
        region = "main";
        D = D1+(x+L0/2)*(D2-D1)/L0;
    end

    % Brownian step
    dX = sqrt(2*D*dt)*randn(1,3);
    if x>=-L0/2 && x<=L0/2 % drift if in bulk
        x_new = x + alpha*((D2-D1)/L0)*dt + dX(1);
    else
        x_new = x + dX(1);
    end
    y_new = y + dX(2);
    z_new = z + dX(3);

    r_new = sqrt(y_new^2 + z_new^2);

    % ---------- Radial reflection ----------
    if region == "main"
        R = R0;
    else
        R = R1;
    end

    y_new = R - abs(mod(r_new + R, 4*R) - 2*R); % reflect if you need to
    z_new=0;
    r_new = y_new;

    % ---------- Axial logic ----------
    % Attempting to leave main cylinder
    if region == "main"
        if x_new < -L0/2
            if r_new <= R1
                % enter left branch
                x = x_new;
            else
                % reflect in x
                x = -L0/2 + (-L0/2 - x_new);
            end
        elseif x_new > L0/2
            if r_new <= R1
                % enter right branch
                x = x_new;
            else
                % reflect in x
                x = L0/2 - (x_new - L0/2);
            end
        else
            x = x_new;
        end

        % In left branch
    elseif region == "left"
        if x_new <= -L0/2 - L1
            T = t;   % absorption
            exitRight = 0;
            output = [T exitRight];
            return
        else
            x = x_new;
        end

        % In right branch
    else
        if x_new >= L0/2 + L2
            T = t;   % absorption
            exitRight = 1;
            output = [T exitRight];
            return
        else
            x = x_new;
        end
    end

    % Update transverse coordinates
    y = y_new;
    z = z_new;

    t = t + dt;
end

T = NaN;
end
