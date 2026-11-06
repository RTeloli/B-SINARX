function [z ,acc , v , displ , u , dz , t] = getSeveralExcitation(varargin)
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
selectedAmplitude = ["15"  "30"];
selectAmp = randi(length(selectedAmplitude),  nExperNeeded, 1);
i = 1;
data = load('Chirp0to40.mat'); 

for l = 1:nExperNeeded 
    [z(:, i),acc(:, i),... 
        v(:, i), displ(:, i), ...
        u(:, i), dz(:, i), t(:, i)] = getCorrrectedData04(data, "n_experiments", 1, "amplitude",selectedAmplitude(selectAmp(l)) );
    i = i+1;

end

z    = z(:); 
acc    = acc(:); 
v    = v(:); 
displ    = displ(:); 
u    = u(:); 
dz    = dz(:); 
t    = t(:); 





end
