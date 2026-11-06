function [z, u, v, x , a, t, params , dZ] = getRandomized(nExperNeeded,varargin)

fc                     = 80;
Fileindex              = 2;
SignalTypes = ["Chirp" ];
Ampls  = ["30","15", "30"  ];
upsample = 0;
% Ampls  = ["15"  ];
for i = 1:length(varargin)
    if isstring(varargin{i})
        switch varargin{i}.lower()
            case "signaltype"
                SignalTypes = varargin{i+1};
            case "upsample"
                upsample = varargin{i+1};


        end
    end
end

indexSignal = randi(length(SignalTypes), nExperNeeded, 1);
indexAmpls  = randi(length(Ampls), nExperNeeded, 1);

z = [];
u = [];
v = [];
x = [];
a = [];
t = [];
dZ = [];
for j = 1:nExperNeeded

    [Z, U, V, X , A, T, params ]   = getDataFolderFiles( ...
        "n_experiments",     1     ,...
        "amplitude",         Ampls(indexAmpls(j)),...
        "cuttoffFreq",       fc,...
        "index",             Fileindex,...
        "signal",            SignalTypes(indexSignal(j)) );
    tselect =  (T<6& T>3);
    Z = Z(tselect);
    U = U(tselect);
    V = V(tselect);
    X = X(tselect);
    A = A(tselect);
    T = T(tselect);

    dz = gradient(Z, T);
    z   =[z; Z];
    u   =[u; U];
    v   =[v; V];
    x   =[x; X];
    a   =[a; A];
    t   =[t  T];
    dZ = [dZ; dz];

end
Ampls(indexAmpls)
h = t(2)-t(1);
fs = 1/h;
fsint = upsample*fs;
hsint = 1/fsint;
tsint = t(1):hsint:t(end);
samplertyp = "cubic";
if logical(upsample)

    z = interp1(t,z , tsint, samplertyp)';
    u = interp1(t,u , tsint, samplertyp)';
    v = interp1(t,v , tsint, samplertyp)';
    x = interp1(t,x , tsint, samplertyp)';
    a = interp1(t,a , tsint, samplertyp)';
    t = tsint;
end


