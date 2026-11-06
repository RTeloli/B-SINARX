function [z ,acc , v , displ , u , dz  ] = getSeveralExcitation02(varargin)
nExperNeeded = 3;
selectedSignal = "chirp";
for i  = 1:length(varargin)
    if isstring(varargin{i})
        switch varargin{i}
            case "n_experiments"
                nExperNeeded = varargin{i+1};
        end
    end
end
selectedAmplitude = ["30" ];
selectAmp = randi(length(selectedAmplitude),  nExperNeeded, 1);
i = 1;
data = load('Chirp0to40.mat'); 

Z     = []; 
Acc   = []; 
V     = []; 
Displ = []; 
U     = []; 
Dz    = []; 
for l = 1:nExperNeeded 
    if ~strcmp(selectedAmplitude(selectAmp(l)), "rand")
    [z, acc, v, displ, u, dz ] = getCorrrectedData04(data, "n_experiments", 1, "amplitude",selectedAmplitude(selectAmp(l)) );
    
    i = i+1;
    else
        [z, acc, v, displ, u, dz ] = getCorrrectedData(data, "n_experiments", 1, "signal",selectedAmplitude(selectAmp(l)) );
        i = i+1;

    end

    Z     = [Z    ; z];
    Acc   = [Acc  ; acc];
    V     = [V    ; v];
    Displ = [Displ; displ];
    U     = [U    ; u];
    Dz    = [Dz   ; dz];

end

z      =  Z    ;
acc    =  Acc  ;
v      =  V    ;
displ  =  Displ;
u      =  U    ;
dz     =  Dz   ; 





end
