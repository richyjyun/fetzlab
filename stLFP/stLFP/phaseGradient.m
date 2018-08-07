function [fx,fy] = phaseGradient(phase,k)
% Summary
%   Calculates the phase gradient of the grid of phases. ASSUMES USE WITH
%   IMAGESC (FLIPPED YDIR)
%
% Inputs
%   phase       Grid of phase values (-pi to pi)
%   k           Distance of neighbors to use for calculation
%
% Outputs
%   fx          Real part of phase gradient (x component. For quiver)
%   fy          Imaginary part of phase gradient (y component)
%
% RJY 07//18/2018

if ~exist('k','var')
    k = 2;
end

gradient = nan(size(phase));

for y = 1:size(phase,1)
    
    for x = 1:size(phase,2)
        
        if(any(isnan(phase(y,x,:)))), continue; end
        
        gradient(y,x,:) = 0;
        
        for i = -k:k
            
            for j = -k:k
                
                if(i==0 && j==0), continue; end
                
                if(x+i<=0 || x+i>size(phase,2) || y+j<=0 || y+j>size(phase,1)), continue; end
                
                if(any(isnan(phase(y+j,x+i,:)))), continue; end
                
                p = squeeze(phase(y+j,x+i,:)-phase(y,x,:));
                
                div = sign(p).*pi.*ones(length(p),1);
                div(div==0) = pi;
                
                p = mod(p,div);
                
                d = sqrt(i^2+j^2);
                
                a = atan2(-j,i); % need negative j just because of how rows/columns are defined
                
                tempGradient(1,1,:) = (p/d)*exp(1j*a);
                
                gradient(y,x,:) = gradient(y,x,:) + tempGradient;
                
            end
            
        end
        
    end
    
end

fx = real(gradient);
fy = -imag(gradient); % need to be negative to make it display properly (imagesc flips yaxis)

end


