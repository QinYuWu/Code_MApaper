function [f] = tlspdf(x,mu,sigma,nu)
f=(gamma((nu+1)./2)./(sigma*sqrt(nu*pi)*gamma(nu/2)))*((nu+((x-mu)./sigma).^2)./nu).^(-(nu+1)./2);
end
