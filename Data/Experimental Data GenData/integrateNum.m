function [x_num, v_num] = integrateNum(a, t, x0, v0)
     
v_num  = cumtrapz(t, a) + v0;  
v_num  = detrend(v_num) ;
x_num  = cumtrapz(t,cumtrapz(t, a)) + x0;
x_num  = detrend(x_num) ;
end