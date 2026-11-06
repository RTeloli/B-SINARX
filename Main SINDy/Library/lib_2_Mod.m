function [solMod ,lmod,   string, mainStringMod] = lib_2_Mod(x, p, vars)

[sol ,l,   string, mainString]                     = lib_2(x, p, vars); 

% [solMod, ia,idup ] = unique(sol', "rows");
[solMod, ia ] = uniquetol(sol','ByRows',true);
 
lmod= l(ia, :); 
mainStringMod = mainString(ia);  
solMod = solMod';   

 
end
