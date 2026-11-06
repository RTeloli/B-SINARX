function [x, v, a,z ,u, h, t, data] = getBoucWendata(nper, options)
arguments
    nper 
    options.type   (1, 1) string {mustBeMember(options.type , ["multisine", "sin","chirp"])} = "multisine"
    options.tspan  (1, :) double = [0 1]
    options.upsamp (1, 1) double = 10
    options.fs     (1, 1) double = 750
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$%%%
% Initial test
%
% Authors: Villani and Teloli - UNESP, Ilha Solteira
% Advisor: Samuel da Silva 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$%%%
 
%% Model Parameters.
extType =  options.type  ;
tspan   =  options.tspan ;
fs      =  options.fs    ;
upsamp  =  options.upsamp;

 
data = struct(); 

data.model.m     = 2;      % mass [kg]
data.model.c     = 10;         % damping coeficient [Ns/m]
data.model.k     = 5e4;        % linear stiffness [N/m]
data.model.alpha = 5e4;        % Bouc-Wen parameter
data.model.beta  = 1e3;        % Bouc-Wen parameter
data.model.gamma = 0.8;        % Bouc-Wen parameter
data.model.delta = -1.1;       % Bouc-Wen parameter
data.model.ni    = 1;          % Bouc-Wen parameter
data.model.wn    = sqrt((data.model.k+data.model.alpha)/data.model.m); % natural frequency
 
%% Time integration parameters.
fmin = 5;               % excitation bandwidth.
fmax = 150;
 
switch extType
    case "sin"
        fmin_selected = fmin + (fmax - fmin)*rand(1,1);
        T = 1/fmin_selected; 
        tspan =[0 T*15];
end
fsint   = fs*upsamp;          % integration sampling frequency.
h       = 1/fsint;            % integration time step.

P = 1;                                % number of excitation periods.
N = tspan(end)/(upsamp*h*(P)) + 1/P;  % number of points per period.
Nint = floor(N*upsamp);               % number of points per period during integration.
t =(0:Nint*P-1)*h;
 

Acond = 50;   
 
switch extType
    case "multisine"
        U = zeros(Nint,1);      % definition of the multisine excitation.
        fres = fsint/Nint;
        exclines = 1:ceil(fmax/fres);
        exclines(exclines < floor(fmin/fres)) = [];


        U(exclines+1) = exp(complex(0,2*pi*rand(size(exclines))));
        u             = 2*real(ifft(U));
        tmp           = Acond*u./std(u);
        u             = repmat(tmp,[P 1]);
    case "sin"
        fmin_selected = fmin + (fmax - fmin)*rand(1,1); 
        u = sin(2*pi*fmin_selected*t);
        tmp           = Acond*u./std(u);
        u             = repmat(tmp,[P 1])';
    case "chirp" 

        u = chirp(t,fmin,t(end)*2,fmax);
        tmp           = Acond*u./std(u);
        u             = repmat(tmp,[P 1])';
    otherwise
        u = zeros(Nint,1); 

end
 
% Time integration.
 
[tmp, tmp0, tmp2, tmp1, ~] = ode45u('bouc', t, [0;0;0], u',data.model);
x  = addnoise(tmp(1:(P)*Nint),nper/100);
v  = addnoise(tmp0(1:(P)*Nint),nper/100);
a  = addnoise(tmp2(1:(P)*Nint),nper/100);
z  = addnoise(tmp1(1:(P)*Nint),nper/100);
u = addnoise(u(1:(P)*Nint),nper/100);

 
end


