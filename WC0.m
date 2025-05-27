function [y] = WC0(k)
%Worst-case approach with mu=1 and sigma=1
for i=1:length(k)
    y(i)=sqrt((k(i)-1)^2/(2*k(i)-1));
end
end

