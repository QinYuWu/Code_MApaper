function [y] = rho_t(k,nu)
%UNTITLED2 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
n=length(k);
for i=1:n
    y(i)=integral(@(x) k.*x.^(k-1).* tinv(x,nu),0,1);
end
end
