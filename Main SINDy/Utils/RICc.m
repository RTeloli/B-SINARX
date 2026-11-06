function [RiCc, k,rss ] = RICc( x_theta, y, beta)
M = size(x_theta, 1);
left = x_theta'*y;
betals =  pinv(x_theta'*x_theta)*left;
RiCc  =zeros(length(beta), size(y, 2));
k =zeros(length(beta), size(y, 2));
n_outputs = size(y, 2);
for i=1:length(beta)

    P = l_0(betals);
    pi_rcc = 2*(log(P) + log(log(P)));
    beta_0 = l_0(beta{i});
    xb = x_theta*beta{i};
    % xb = closedloop(obj, x, beta{i})'; 
    for j = 1:n_outputs
        y_hat =xb(:, j);
        y_ls =x_theta*betals(:,j);
        y_true = y(:, j);
        Rss = norm(y_hat - y_true,2) ; 
        sigma_hat = norm(y_ls-y_true, 2)/(M-P(j));

        RiCc(i, j) = Rss + pi_rcc(j)*sigma_hat*beta_0(j);
        k(i, j)    = beta_0(j);
        rss(i, j) = Rss;

    end
end

% for i = 1:size(y, 2)
%     figure(i*10)
%     semilogy(k(:, i), abs(RiCc(:, i)))
%     grid on
% 
%     xlabel("k")  
%     ylabel("RICc")
% end
end
function l0 = l_0(beta)
l0 = sum(abs(beta)>0);
end