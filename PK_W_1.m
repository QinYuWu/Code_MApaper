function [y] = PK_W_1(p,epsilon,k,n)   %p:moment; epsilon:uncertainty; k:order of PD
y=0;
alpha=0:(1/n):1;
for i=1:(n-1)
   y=y+qW_1(alpha(i+1),p,epsilon)*(alpha(i+1)^k-alpha(i)^k);
end
y=y+qW_1(0.9999,p,epsilon)*(1-(1-1/n)^k);
end



