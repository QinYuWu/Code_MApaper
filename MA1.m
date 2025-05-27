function [y] = MA1(k)
%Worst-case approach with mu=1 and sigma=1
for i=1:length(k)
    y(i)=sqrt(pi)*gamma(k(i)+1/2)/gamma(k(i));
end
end
