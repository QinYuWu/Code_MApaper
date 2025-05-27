function [y] = object(x,k,alpha)
if x<=norminv(alpha,0,1)
   y=0;
else
beta=normcdf(x,0,1);
y=quadl(@(s) (x-norminv(s,0,1)).^k,alpha,beta);
end
