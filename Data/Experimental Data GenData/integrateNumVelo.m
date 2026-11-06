function x_num = integrateNumVelo(v, t, x0)
     
 
x_num  = cumtrapz(t,v) + x0;
x_num  = detrend(x_num) ;
end