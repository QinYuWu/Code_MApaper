function [q] = qW_1(alpha,k,epsilon)
n=length(alpha);
for i=1:n
    q(i)=fsolve(@(x) object(x,k,alpha(i))-epsilon^k,norminv(alpha(i),0,1)+1);
end
end


