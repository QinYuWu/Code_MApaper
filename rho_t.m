function [y] = rho_t(k,nu)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
n=length(k);
for i=1:n
    y(i)=integral(@(x) k.*x.^(k-1).* tinv(x,nu),0,1);
end
end
