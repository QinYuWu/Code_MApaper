function [q] = qt_sp1(data,u)
%UNTITLED27 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
for i=1:length(u)
    q(i)=max([fsolve(@(x) normcdf(x,-0.0023,0.0225)-u(i),0),fsolve(@(x) tlscdf(x,-0.00253866,0.0141682,3.08178)-u(i),0),qt_cdf(data,u(i))]);
end
end


