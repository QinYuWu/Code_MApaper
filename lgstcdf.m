function [F] = lgstcdf(x,mu,sigma)
%UNTITLED4 此处显示有关此函数的摘要
%   此处显示详细说明
for i=1:length(x)
F(i)=quadl(@(y) lgstpdf(y,mu,sigma),-1000,x(i));
end
end