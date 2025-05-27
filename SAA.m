function [y] = SAA(data,w,k)
%PW_k with nominal data
loss=data*w';
loss=sort(loss);
n=length(loss);
a=0:(1/n):1;
for i=1:length(k)
for j=1:n
    delta(j)=a(j+1)^k(i)-a(j)^k(i);
end
y(i)=sum(loss'.*delta);
end
end

