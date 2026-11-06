function [X_train, y_train, testData, inputvariables, outputvariables, delays, arecur, vars, params] = getDataForSINDy3()


extType = "multisine";
z = []; 
u = []; 
v = []; 
x = []; 
a = []; 
t = []; 
noiselvl = 0.0;
for i = 1:5
    [X  ,V , Z, U, T ] = getBWModelData(extType);
    z = [z Z ]; 
    u = [u U ]; 
    v = [v V ]; 
    x = [x X ];   
end
 
Y        = [x;  v ; z;  u  ]';
vars     = {"x", "v", "z", "u"} ;
delays   = [1 1 1 1] ;
arecur   = logical([1 1 1  0 ]);
[input, output, ...
    inputvariables, outputvariables] = getTrainingData3(Y, ...
    "vars", vars, ...
    "recurIndex", arecur, ...
    "delays", delays, ...
    "noise", noiselvl);


nExperNeeded           = 1;
excitationLevel        = "30";
Fileindex = 2;
SignalType             = "Chirp" ;

extType = "chirp";
[xtest  ,vtest , ztest, utest, t, params  ] = getBWModelData(extType);
 
Ytest = [xtest; vtest; ztest; utest ]';

[input_test, output_test]  = getTrainingData3(Ytest, ...
    "vars"       , vars, ...
    "recurIndex" , arecur, ...
    "delays"     , delays);




X_train        = input' ;
y_train        = output';
%
testData.xtest = input_test';
testData.ytest = output_test';

end

function [xode ,vode, zode, force, t, params] = getBWModelData(extType)

fmin   = 5;               % excitation bandwidth.
fmax   = 40;
fs     = 1024  ; 
upsamp = 1; 
tspan  = [0,10]; 

switch extType
    case "sin"
        fmin_selected = fmin + (fmax - fmin)*rand(1,1);
        T = 1/fmin_selected;
        tspan =[0 T*15];
end
fsint    = fs*upsamp;          % integration sampling frequency.
h        = 1/fsint;            % integration time step.

params.h = h;     
P        = 1;
N        = tspan(end)/(upsamp*h*(P)) + 1/P;              
Nint     = floor(N*upsamp);         
t        = (0:Nint*P-1)*h; 
Acond    = 200*(2-  rand(1, 1));         
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

        u = chirp(t,fmin,t(end) ,fmax);
        tmp           = Acond*u./std(u);
        u             = repmat(tmp,[P 1])';
    otherwise
        u = zeros(Nint,1);

end
delta_vol = 341;
alpha_vol = 5.000292979243245e+02;
k2 = 0;
Ac = 49.5;
omega_n = 18.80*2*pi;
m = 1-alpha_vol/omega_n^2;
zeta = 0.00446;     %
gamma_vol = 330.01;


params.m = 1; 
params.c = omega_n*2*zeta; 
params.k = omega_n^2; 
force = u';
tforce = t;
x0 = [0; 0];

[xode ,vode , ~, zode, ~] = ode45u('bouc2', tforce, ...
    [x0;0], ...
    force, ...
    [omega_n ...
    zeta ...
    k2 ...
    alpha_vol ...
    1 ...
    gamma_vol ...
    delta_vol ...
    1 ...
    (omega_n^2-alpha_vol)]);
end