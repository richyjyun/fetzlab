%% Helper function for performing granger causality in pre/cond/post. 
% For cs_CohGC2.

function power = doGC(params,trig1,trig2,trig3,inds,Filter,chnsL,chnsR,chnm,fig,line)
        power = []; % (epoch,beta vs gamma, rchn,lchn,l2r vs r2l)

        for t = 1:3
            X = [];
            
            switch t
                case 1
                    trig = trig1;
                case 2
                    trig = trig2;
                case 3
                    trig = trig3;
            end     
            
            trialinds = repmat(trig', length(inds), 1) + repmat(inds(:), 1, size(trig,1));
            trialinds(:,floor(trialinds(1,:))<=0) = []; trialinds(:,floor(trialinds(end,:))>length(Filter)) = [];
            
            params.ntrials   = length(trig);     % number of trials
            params.nobs      = length(inds);   % number of observations per trial
            
            for i = 1:length(chnsR)
                d = Filter(:,chnsR(i));%filtfilt(b,a,double(data(:,chnsR(i)))); %comb notch
                Snips = u.meanSubtract(d(floor(trialinds)),params);
                X(1,:,:) = Snips;
                for j = 1:length(chnsL)
                    d2 = Filter(:,chnsL(j));%filtfilt(b,a,double(data(:,chnsL(j)))); %comb notch
                    Snips2 = u.meanSubtract(d2(floor(trialinds)),params);
                    X(2,:,:) = Snips2;
                    
                    disp([num2str(i),'-',num2str(j)]);
                    
                    [gc,f] = a.GC(X,params); % custom helper function for using MVGC toolbox
                    
                    power(t,1,i,j,1) = sum(gc(2,f>=params.beta(1)&f<=params.beta(2)));
                    power(t,1,i,j,2) = sum(gc(1,f>=params.beta(1)&f<=params.beta(2)));
                    power(t,2,i,j,1) = sum(gc(2,f>=params.gamma(1)&f<=params.gamma(2)));
                    power(t,2,i,j,2) = sum(gc(1,f>=params.gamma(1)&f<=params.gamma(2)));                
                end
            end
        end
        
        set(0, 'CurrentFigure', fig);
        for i = 1:4*length(chnsR)*length(chnsL)
            subplot(length(chnsL),4*length(chnsR),i)
            if(mod(i-1,4*length(chnsR)) < 2*length(chnsR))
                freq = 1;
            else
                freq = 2;
            end
            
            if(mod(i-1,2*length(chnsR))< length(chnsR))
                dir = 1;
            else
                dir = 2;
            end
            
            r = mod(i-1,length(chnsR))+1;
            l = floor((i-1)/(4*length(chnsR))) + 1;
            
            hold on; plot(power(:,freq,r,l,dir),line); hold off;
            if r == 1 && freq ==1 && dir == 1
                ylabel(chnm(chnsL(l)));
            end
            if l == 1
                title(chnm(chnsR(r)));
            end
        end
end