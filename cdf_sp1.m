function [F] = cdf_sp1(data,x)
for i=1:length(x)
    F(i)=min([normcdf(x(i),-0.0023,0.0225),tlscdf(x(i),-0.00253866,0.0141682,3.08178),cdf(data,x(i))]);
end
end

