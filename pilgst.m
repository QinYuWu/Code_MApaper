function [y] = pilgst(x,mu,sigma)
%UNTITLED9 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
for i=1:length(x)
y(i)=quadl(@(z) 1-lgstcdf(z,mu,sigma),x(i),0.15);
end
end

