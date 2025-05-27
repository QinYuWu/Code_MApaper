function [y] = WC_pi(data,x)
for i=1:length(x)
    y(i)=max([picdf(x(i),data),quadl(@(z) 1-normcdf(z,mean(data),sqrt(var(data))),x(i),1),pitls(x(i),-0.00253866,0.0141682,3.08178),pilgst(x(i),-0.00248069,0.011185)]); 
end
end
