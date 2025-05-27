function [y] = ES(alpha,F)
%UNTITLED18 此处显示有关此函数的摘要
%   此处显示详细说明
x_l=fsolve(@(z) F(z)-alpha,0.035);

y=1/(1-alpha)*(1-x_l*alpha-quadl(@(z) F(z),x_l,1));

end
