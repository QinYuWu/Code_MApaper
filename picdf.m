function [y] = picdf(x,data)
%UNTITLED12 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
n=length(data);
data=sort(data);
for i=1:length(x)
    for j=1:n
        if x(i)>=data(j)
            SL(j)=0;
        else SL(j)=data(j)-x(i);
        end
        y(i)=mean(SL);
    end
end



