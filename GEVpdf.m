function [f] = GEVpdf(x,mu,sigma,k)
f=(1./sigma).*exp(-(1+k*(x-mu)./sigma).^(-1./k)).*(1+k*(x-mu)./sigma).^(-1-1./k);
end

