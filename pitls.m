function [y] = pitls(x,mu,sigma,nu)
%UNTITLED9 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
for i=1:length(x)
y(i)=quadl(@(z) 1-tlscdf(z,mu,sigma,nu),x(i),0.5);
end
end

