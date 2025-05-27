function [y] = MA0(k)
%Model aggregation approach with mu=1 and sigma=1
g=@(x,s) x.*(1-x).^(-0.5).*(1+x).^(s-1.5);
for i=1:length(k)
y(i)= k(i).*0.5.^k(i)*integral(@(x)g(x,k(i)),-1,1);
end
end

