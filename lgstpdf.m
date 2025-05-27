function [f] = lgstpdf(x,mu,sigma)
f=exp((x-mu)./sigma)./(sigma*(1+exp((x-mu)./sigma)).^2);
end

