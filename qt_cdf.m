function [q] = qt_cdf(data,u)
n=length(data);
data=sort(data);
for i=1:length(u)
    j=1;
       while(u(i)>j/n)
           j=j+1;
       end
q(i)=data(j);
end
end

