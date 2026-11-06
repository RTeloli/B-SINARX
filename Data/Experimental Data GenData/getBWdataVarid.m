function [yf, dyf, ddyf, zf,u,h, t] = getBWdataVarid(nper, varargin)
% 
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$%%%
% Initial test
%
% Authors: Villani and Teloli - UNESP, Ilha Solteira
% Advisor: Samuel da Silva 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$%%%
% clear all; close all; clc
%% Model Parameters.
extType = "multisine";
tspan   = [0 1];
fs      = 750;               % working sampling frequency.
upsamp  = 10;            % upsampling factor to ensure integration accuracy.
dynamicsType = "bouc";
for i = 1:length(varargin)
    if isstring(varargin{i})
        switch varargin{i}
            case "type"
               extType = varargin{i+1};
            case "tspan"
               tspan   = varargin{i+1};
            case "upsamp"
               upsamp  = varargin{i+1};
            case "dynaType"
                dynamicsType = varargin{i+1};
        end
    end
end

data = struct();
model = struct();

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
fsint   = fs*upsamp;      % integration sampling frequency.
h       = 1/fsint;            % integration time step.

P = 1;                  % number of excitation periods.
N = tspan(end)/(upsamp*h*(P)) + 1/P;               % number of points per period.
Nint = floor(N*upsamp);        % number of points per period during integration.
t =(0:Nint*P-1)*h;

% fmin = 35.58;               % excitation bandwidth.
% fmax = 35.58;

Acond = 50; 10:10:80;                 % excitation amplitude.

Pfilter = 1;            % extra period to avoid edge effects during low-pass filtering (see line 59).
% P = P + Pfilter;
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
        % fmin = 200;
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
% 
% Low-pass filtering and downsampling.
drate = factor(upsamp);        % prime factor decomposition.
y = []; z= []; yf = []; uf = []; zf = [];
% nper = 1; % dB noise

switch dynamicsType
    case "bouc"

        [tmp, tmp2, tmp3, tmp1, ~] = ode45u('bouc', t, [0;0;0], u',data.model);
        yf   = addnoise(tmp(1:(P)*Nint),nper/100);
        dyf  = addnoise(tmp2(1:(P)*Nint),nper/100);
        ddyf = addnoise(tmp3(1:(P)*Nint),nper/100);
        zf   = addnoise(tmp1(1:(P)*Nint),nper/100);
        u    = addnoise(u(1:(P)*Nint),nper/100);  
 
    case "duff"
        m     = 1; 
        c     = 10; 
        k     = 100; 
        gamma = 0.2; %

        data.model.gamma = gamma;    
         data.model.m    = m ;
         data.model.c    = c ;
         data.model.k    = k ;  
        [y, dy, ddy, z,~] = getDataDuffing(u', t);%%
        yf    = addnoise(y(1:(P)*Nint),nper/100);
        dyf   = addnoise(dy(1:(P)*Nint),nper/100);
        ddyf  = addnoise(ddy(1:(P)*Nint),nper/100);
        zf    = addnoise(z(1:(P)*Nint),nper/100);
            
        u   = addnoise(u(1:(P)*Nint),nper/100);
 
end
%% Low-pass filtering and downsampling.

% drate = factor(upsamp);        % prime factor decomposition.
% 
% for k=1:length(drate)
%     yf = decimate(yf,drate(k),'fir'); 
%     zf = decimate(zf,drate(k),'fir');    
% end %k
% 
% u = downsample(u,upsamp);

%% Removal of the last simulated period to eliminate the edge effects due to the low-pass filter.
% 
% yf = yf(1:(P-1)*N);
% zf = zf(1:(P-1)*N);
% u = u(1:(P-1)*N,:);
% P = P-1;

% z = u' - (ddy*2 + 10*dy + 5e4*yf);

end


