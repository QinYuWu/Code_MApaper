function [F] = tlscdf(x,mu,sigma,nu)
%UNTITLED4 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
for i=1:length(x)
F(i)=quadl(@(y) tlspdf(y,mu,sigma,nu),-1000,x(i));
end
end

