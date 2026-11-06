function [theta1_reduced, theta2_reduced] = reduceTheta(theta1, theta2, xi_0)
premRemove          = all(xi_0'==0);
theta1_reduced = theta1(:, ~premRemove);    
theta2_reduced = theta2(:, ~premRemove);  
 

end