function [dxdt, fi] = bouc2(t, x, u, p)

% state vector
v_disp = x(1);
v_vel  = x(2);
v_z    = x(3);

% parameters
omega_n = p(1);      % mass [kg]
zeta    = p(2);         % damping coeficient [Ns/m]
k2      = p(3);           % linear stiffness [N/m]
alpha   = p(4);        % Bouc-Wen parameter
beta    = p(5);         % Bouc-Wen parameter
gamma   = p(6);        % Bouc-Wen parameter
delta   = p(7);        % Bouc-Wen parameter
ni      = p(8);           % Bouc-Wen parameter
k       = p(9); 

fi = v_z;

% state equations
dxdt       = [v_vel;
             (u  - v_vel*2*zeta*omega_n - v_disp*k- v_z - k2*v_disp^2);
              alpha*v_vel - beta*(gamma*abs(v_vel)*abs(v_z)^(ni-1)*v_z+delta*v_vel*abs(v_z)^ni);
              ];
end

