function [F] = lgstcdf(x,mu,sigma)
%UNTITLED4 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
for i=1:length(x)
F(i)=quadl(@(y) lgstpdf(y,mu,sigma),-1000,x(i));
end
end