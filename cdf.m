function [F] = cdf(data,x)
%UNTITLED20 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
n=length(data);
data=sort(data);
for i=1:length(x)
if x(i)<=data(1)
    F(i)=0;
elseif x(i)>=data(n)
    F(i)=1;
else 
    j=1;
    while(x(i)>data(j))
        j=j+1;
    end
    F(i)=j/n;
end
end
