function [y] = rho_norm(k)
%UNTITLED2 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
n=length(k);
for i=1:n
    y(i)=integral(@(x) k.*x.^(k-1).* norminv(x,0,1),0,1);
end
end

