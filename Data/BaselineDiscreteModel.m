function [xi_s, expressions, vars] = BaselineDiscreteModel(options)
arguments
    options.numerical    (1,1) string {mustBeMember(options.numerical,    ["boucwen0","linearbw", "boucwen1", "boucwen2"])}
    options.experimental (1,1) string {mustBeMember(options.experimental, ["linear_no", "linear_ex"])}
    options.parameters   struct = struct()
    options.timestep     (1,1) double {mustBePositive} = 1/1024
end

% ---- Exactly one data mode must be given ----
hasNum = isfield(options, "numerical");
hasExp = isfield(options, "experimental");

if hasNum && hasExp
    error('BaselineDiscreteModel:bothModes', ...
        ['Select only one data mode.\nMake sure you do:\n' ...
         '   1. BaselineDiscreteModel(numerical="boucwen1")  or\n' ...
         '   2. BaselineDiscreteModel(experimental="linear_ex")']);
elseif hasNum
    datamode   = "numerical";
    modeselect = options.numerical;
elseif hasExp
    datamode   = "experimental";
    modeselect = options.experimental;
else
    datamode   = "experimental";    % default
    modeselect = "linear_no";
end

dt = options.timestep;

switch datamode
    case "numerical"
        switch modeselect
            case "boucwen0"
                % Bouc-Wen Linear Baseline
                m     = 2;          % mass [kg]
                c     = 10;         % damping coeficient [Ns/m]
                k     = 5e4;        % linear stiffness [N/m]


                xi_s = zeros(3, 3);

                xi_s(1, 1) = 0;
                xi_s(2, 1) = 0;

                xi_s(1, 2) = 0;                  % x_k
                xi_s(2, 2) = 0;               % v_k
                xi_s(3, 2) = 0;

                expressions = {};
                expressions{1} =  "x1";
                expressions{2} =  "v1";
                expressions{3} =  "u1";

                vars = {"x1", "v1", "z1", "u1"};
            case "linearbw"
                % Bouc-Wen Linear Baseline
                m     = 2;          % mass [kg]
                c     = 10;         % damping coeficient [Ns/m]
                k     = 5e4;        % linear stiffness [N/m]


                xi_s = zeros(3, 3);

                xi_s(1, 1) = 1;
                xi_s(2, 1) = dt;

                xi_s(1, 2) = - dt*k/m;                  % x_k
                xi_s(2, 2) =  1 - dt*c/m;               % v_k
                xi_s(3, 2) = -1*dt/m;

                expressions = {};
                expressions{1} =  "x1";
                expressions{2} =  "v1";
                expressions{3} =  "u1";

                vars = {"x1", "v1", "z1", "u1"};
            case "boucwen1"
                % Bouc-Wen Partial Nonlinear Knowledge |v|z
                m     = 0.2;      % mass [kg]
                c     = 10;         % damping coeficient [Ns/m]
                k     = 5e4;        % linear stiffness [N/m]
                alpha = 5e4;        % Bouc-Wen parameter
                beta  = 1e3;        % Bouc-Wen parameter
                gamma = 0.8;        % Bouc-Wen parameter


                xi_s = zeros(5, 3);

                xi_s(1, 1) = 1;
                xi_s(2, 1) = dt;

                xi_s(1, 2) = - dt*k/m;                  % x_k
                xi_s(2, 2) =  1 - dt*c/m;               % v_k
                xi_s(3, 2) = -1*dt/m;                   % z_k
                xi_s(4, 2) =  -dt/m;                    % u_k

                xi_s(2, 3) = dt*alpha;
                xi_s(3, 3) = 1;
                xi_s(5, 3) = -dt*gamma*beta;

                expressions = {};
                expressions{1} =  "x1";
                expressions{2} =  "v1";
                expressions{3} =  "z1";
                expressions{4} =  "u1";
                expressions{5} =  "abs(v1)z1";

                vars = {"x1", "v1", "z1", "u1"};


            case "boucwen2"
                % Bouc-Wen Full Nonlinear Knowledge
                m     = 2;          % mass [kg]
                c     = 10;         % damping coeficient [Ns/m]
                k     = 5e4;        % linear stiffness [N/m]
                alpha = 5e4;        % Bouc-Wen parameter
                beta  = 1e3;        % Bouc-Wen parameter
                gamma = 0.8;        % Bouc-Wen parameter
                delta = -1.1;       % Bouc-Wen parameter



                xi_s = zeros(6, 3);

                xi_s(1, 1) = 1;
                xi_s(2, 1) = dt;

                xi_s(1, 2) = - dt*k/m;                  % x_k
                xi_s(2, 2) =  1 - dt*c/m;               % v_k
                xi_s(3, 2) = -1*dt/m;                   % z_k
                xi_s(4, 2) =   dt/m;                    % u_k

                xi_s(2, 3) = dt*alpha;
                xi_s(3, 3) = 1;
                xi_s(5, 3) = -dt*gamma*beta;
                xi_s(6, 3) = -dt*delta*beta;

                expressions = {};
                expressions{1} =  "x1";
                expressions{2} =  "v1";
                expressions{3} =  "z1";
                expressions{4} =  "u1";
                expressions{5} =  "|v1|z1";
                expressions{6} =  "v1|z1|";

                vars = {"x1", "v1", "z1", "u1"};


        end
    case "experimental" 
        switch modeselect
            case "linear_ex"
                % Bolted Joint Linear Baseline
                omega_n = 18.80*2*pi;
                eta     = 0.446/100;
                k       = omega_n^2;           
                xi_s = zeros(4, 3);
                xi_s(1, 1) = 1;
                xi_s(2, 1) = dt;
                xi_s(1, 2) = -dt*k;                   % x_k
                xi_s(2, 2) = 1 - 2*dt*eta*omega_n;    % v_k    
                xi_s(3, 2) = -dt;                     % z_k
                xi_s(4, 2) =  dt;                     % u_k

                expressions = {"x1", "v1", "z1", "u1"};
                vars        = {"x1", "v1", "z1", "u1"};

            case "linear_no"
                % Bolted Joint No Baseline knowledge
                xi_s        = zeros(4, 3);
                expressions = {"x1", "v1", "z1", "u1"};
                vars        = {"x1", "v1", "z1", "u1"};
        end
 
end

end

 