function [fx,fy] = phaseGradient(phase,k)

    if ~exist('k','var')
        k = 2;
    end

    gradient = nan(size(phase));
    
    for y = 1:size(phase,1)
        
        for x = 1:size(phase,2)
            
            if(isnan(phase(y,x))), continue; end
            
            gradient(y,x) = 0;
            
            for i = -k:k
                              
                for j = -k:k
                    
                    if(i==0 && j==0), continue; end
                    
                    if(x+i<=0 || x+i>size(phase,2) || y+j<=0 || y+j>size(phase,1)), continue; end
                    
                    if(isnan(phase(y+j,x+i))), continue; end
                    
                    p = phase(y+j,x+i)-phase(y,x);
                    
                    if(sign(p))
                        div = -pi;
                    else
                        div = pi;
                    end
                    
                    p = mod(p,div);
                    
                    d = sqrt(i^2+j^2);
                    
                    a = atan2(-j,i); % need negative j just because of how rows/columns are defined 

                    gradient(y,x) = gradient(y,x) + (p/d)*exp(1j*a);
                                        
                end
                
            end
                        
        end
        
    end
    
    fx = real(gradient);
    fy = -imag(gradient); % need to be negative to make it display properly (imagesc flips yaxis)
    
end


