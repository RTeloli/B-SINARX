
function [v  u, t, h] = getExperData04(data,varargin )
nExperNeeded = 2; 
amplitude = "30"; 
 for i  = 1:length(varargin)
    if isstring(varargin{i})
        switch varargin{i}
            case "n_experiments"
                nExperNeeded = varargin{i+1};
            case "amplitude"
                amplitude = varargin{i+1};

        end
    end
end

% data = load('Chirp0to40.mat'); 
dataFieldNames =fieldnames(data); 
dataSelected = dataFieldNames(contains(dataFieldNames, amplitude)); 
data =data.(string(dataSelected)); 

% Pre processing
u_full   = data.force;
v_full   = data.laser; 
%%
Nruns = 52; 
Ntimestep = length(v_full)/Nruns; 
fs = 1024; 
h = 1/fs; 
t = 0:h:(Ntimestep-1)*h; 
u   =  reshape(u_full,Ntimestep,Nruns);    %uc_chirp0p20 = uc_chirp0p20(:,1);
v   =  reshape(v_full,Ntimestep,Nruns);    %vel_chirp0p20 = vel_chirp0p050(:,1);    
 
selectedIdnex = randperm(Nruns,nExperNeeded); 
u= u(:,selectedIdnex ); v= v(:,selectedIdnex );  
  


end
 