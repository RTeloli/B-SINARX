function [solMod ,lmod,   string, mainStringMod] = lib6Mod(x, p, vars)

[sol ,l,   string, mainString]                     = lib6(x, p, vars); 

% [solMod, ia,idup ] = unique(sol', "rows");
[solMod, ia ] = uniquetol(sol','ByRows',true);
 
lmod= l(ia, :); 
mainStringMod = mainString(ia);  
solMod = solMod';   

 
end
