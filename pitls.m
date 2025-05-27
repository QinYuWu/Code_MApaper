function [y] = pitls(x,mu,sigma,nu)
%UNTITLED9 此处显示有关此函数的摘要
%   此处显示详细说明
for i=1:length(x)
y(i)=quadl(@(z) 1-tlscdf(z,mu,sigma,nu),x(i),0.5);
end
end

