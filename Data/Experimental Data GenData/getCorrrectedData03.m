function [z,acc, v, displ, u, dz, t] = getCorrrectedData03(varargin) 
nExperNeeded = 1;  
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
[ v, u, t, dt] = getExperData03("n_experiments", nExperNeeded, "amplitude", amplitude);
  
Ac =  49.5;

fc     = 90; 
tselec = t<800;   
v      = v(tselec, : );  
u      = Ac*u(tselec, :); 
t      = t(tselec); 
displ   = integrateNumVelo(v,t, 0);

acc  = zeros(size(v));
for i = 1:size(v, 2)
    acc (:, i)   = gradient(v(:, i),t );
end
% acc = acc; 
mu_xi   =  0.444 ;
co_xi   = 0.57/100;
mu_fn   =  18.80; 
wn = 2*pi*mu_fn;
alpha = 5.000292979243245e+02;
k = wn^2 - alpha; 
m = 1-alpha/wn^2; 

displ = butterWorthFilter(displ, dt, fc); 
v     = butterWorthFilter(v   , dt, fc); 
acc   = butterWorthFilter(acc , dt, fc); 
u     = butterWorthFilter(u   , dt, fc); 

z  = u   -  (acc + mu_xi*2*wn*v +wn^2*displ)  ;



% z  = u -  acc ;
% z  = butterWorthFilter(z, dt, fc); 
dz = gradient(z')/dt; 
dz = dz'; 
%%
% upSamRatio = 1; 
% v    = upsampler(v   , t, upSamRatio); 
% acc  = upsampler(acc , t, upSamRatio); 
% u    = upsampler(u   , t, upSamRatio); 
% z    = upsampler(z   , t, upSamRatio); 
% dz   = upsampler(dz  , t, upSamRatio); 
 
% dz   = butterWorthFilter(dz  , dt, fc); 
% [disp, t]= upsampler(disp, t, upSamRatio); 

displ = displ(:);
v    = v(:);
acc  = acc(:);
u    = u(:);
z    = z(:);
dz   = dz(:); 

end