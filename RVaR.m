function [y] = RVaR(alpha,beta,F)
%UNTITLED18 此处显示有关此函数的摘要
%   此处显示详细说明
x_l=fsolve(@(z) F(z)-alpha,0.035);
for i=1:length(beta)
x_u(i)=fsolve(@(z) F(z)-beta(i),0.035);
end
for i=1:length(beta)
y(i)=1/(beta(i)-alpha)*(x_u(i)*beta(i)-x_l*alpha-quadl(@(z) F(z),x_l,x_u(i)));
end
end

