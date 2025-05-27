function [y] = pilgst(x,mu,sigma)
%UNTITLED9 此处显示有关此函数的摘要
%   此处显示详细说明
for i=1:length(x)
y(i)=quadl(@(z) 1-lgstcdf(z,mu,sigma),x(i),0.15);
end
end

