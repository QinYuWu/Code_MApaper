function x = sroot(f,a,b)        %solve root by Newton's method
if f(a)==0
    x=a;
else  if f(b)==0
       x=b;
    else if f(a)>0
            
            while(abs(b-a)>=1e-6)
                c=(a+b)/2; 
                if f(c)>0
                    a=c;
                else if f(c)<0
                      b=c;
                    else 
                        x=c;
                 break;
                    end
                end
            end
            
      else 
            
            while(abs(b-a)>=1e-6)
                c=(a+b)/2; 
                if f(c)>0
                    b=c;
                else if f(c)<0
                      a=c;
                    else 
                        x=c;
                 break;
                     end
                end
            end
            
          end
       end
end
x=c;
end


