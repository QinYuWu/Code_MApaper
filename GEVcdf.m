function [F] = GEVcdf(x,mu,sigma,k)
for i=1:length(x)
F(i)=quadl(@(y) GEVpdf(y,mu,sigma,k),-100,x(i));
end
end

