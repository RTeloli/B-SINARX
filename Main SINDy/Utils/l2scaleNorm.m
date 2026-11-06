function [theta_thilda, H]      = l2scaleNorm(theta)
p = size(theta, 2);
for i = 1:p

    H(i, i) = sqrt(theta(:, i)'*theta(:, i));
end    
% H = eye(size(H));
theta_thilda = H\theta';  
end
