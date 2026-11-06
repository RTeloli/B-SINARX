 function y                            = predict_CL(obj, x, xi)
        theta_X = obj.libfunction(x', obj.maxpower);
        theta_X = obj.scaleMat\theta_X';
        theta_X = theta_X'; 
        y = theta_X*xi;
    end
