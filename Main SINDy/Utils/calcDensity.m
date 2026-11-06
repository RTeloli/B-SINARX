function density= calcDensity(X)
allZeros = sum((X==0));
allsize =  sum((X~=0))+sum((X==0));

spar = allZeros./allsize;
density = 1-spar;
end
