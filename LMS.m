
function [yn,W,en]=LMS(xn,dn,M,mu,itr)

en = zeros(itr,1);             % error
W  = zeros(M,itr);            
% iteration
for k = M:itr                  % kth iteration
    x = xn(k:-1:k-M+1);         
    y = W(:,k-1).' * x;         
    en(k) = dn(k) - y ;        % error of kth iteration
     
    W(:,k) = W(:,k-1) + 2*mu*en(k)*x;
end
yn = inf * ones(size(xn));
for k = M:length(xn)
    x = xn(k:-1:k-M+1);
    yn(k) = W(:,end).'* x;
end
