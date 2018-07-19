% function PlotLFPs(SLNeuro,SL,fname)

fname = fullfile('F:\Dropbox\repos\abogaard\efetz\U\packets', 'LFPsKato.ps');
if(exist(fname))
    delete(fname);
end

dates = extractfield(SL,'Date'); window = 0.5;

for i = 1:length(SLNeuro)
    
    f = []; f(1) = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]); set(f(1),'visible','off');
    
    % find proper SL, get accel snip for left and right (pre)
    SLInd = find(strcmp(dates,SLNeuro(i).Date));
    trig = SL(SLInd).lefttrials+SL(SLInd).rts_l; trig = trig(SL(SLInd).lefttrials(:,2)<SL(SLInd).trig1(1));
    [LSnip,accelinds] = getAccelSnips(SL(SLInd).accel_raw_l,trig,SL(SLInd).fs,window);
    trig = SL(SLInd).righttrials+SL(SLInd).rts_r; trig = trig(SL(SLInd).righttrials(:,2)<SL(SLInd).trig1(1));
    [RSnip,~] = getAccelSnips(SL(SLInd).accel_raw_r,trig,SL(SLInd).fs,window);
    
    inds = SLNeuro(i).tneuro > -window & SLNeuro(i).tneuro < window;
    
    % loop through channels
    for chn = 1:size(SLNeuro(i).Lbeta,2)
        fInd = chn*2-1;
        if(fInd > 36)
            if(length(f)<2)
                f(2) = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]); set(f(2),'visible','off');
            end
            fInd = chn*2-1-36;
        end
                        
        % loop through epochs
        for epoch = 1:3
            if(epoch <= size(SLNeuro(i).Lbeta,1))
                subaxis(6, 6, fInd, 'padding', .01, 'spacing', .01, 'spacingvert', .03)
                hold on; yyaxis left; plot(SLNeuro(i).tneuro(inds),squeeze(SLNeuro(i).Lbeta(epoch,chn,inds)),'linewidth',2);
                hold on; yyaxis right; plot(SLNeuro(i).tneuro(inds),squeeze(SLNeuro(i).Lgamma(epoch,chn,inds)),'linewidth',2);
                if(epoch == size(SLNeuro(i).Lbeta,1))
                    yl = ylim; hold on; plot(accelinds,(LSnip/max(LSnip))*(yl(2)-yl(1))/4+yl(1),'k-','linewidth',2);
                    ylim(yl);
                end
                title([SLNeuro(i).chnm(SLNeuro(i).chID(chn)),' Left'],'fontsize',7)
                set(gca,'fontsize',7); xlim([-window,window])
            end
            
            if(epoch <= size(SLNeuro(i).Rbeta,1))
                subaxis(6, 6, fInd+1, 'padding', .01, 'spacing', .01, 'spacingvert', .03)
                hold on; yyaxis left; plot(SLNeuro(i).tneuro(inds),squeeze(SLNeuro(i).Rbeta(epoch,chn,inds)),'linewidth',2);
                hold on; yyaxis right; plot(SLNeuro(i).tneuro(inds),squeeze(SLNeuro(i).Rgamma(epoch,chn,inds)),'linewidth',2);
                if(epoch == size(SLNeuro(i).Lbeta,1))
                    yl = ylim; hold on; plot(accelinds,(RSnip/max(RSnip))*(yl(2)-yl(1))/4+yl(1),'k-','linewidth',2);
                    ylim(yl);
                end
                title([SLNeuro(i).chnm(SLNeuro(i).chID(chn)),' Right'],'fontsize',7)
                set(gca,'fontsize',7); xlim([-window,window])
            end
        end
    end
    
    for t = 1:length(f)
        set(0, 'CurrentFigure', f(t))
        a = axes; t1 = title([SLNeuro(i).Date,'; Left=Beta, Right=Gamma']);
        a.Visible = 'off'; t1.Visible = 'on';
        print(f(t), '-dpsc2', fname, '-append')
        close(f(t))
    end
end

function [Snip,inds] = getAccelSnips(Accel,trig,fs,window)

Accel = u.FilterAccForMovementTimes(Accel,fs, 'richardson');
inds = (-window*fs:1:window*fs);
trig(isnan(trig)) = []; trig = trig*fs/1000;
trialinds = repmat(trig', length(inds), 1) + repmat(inds(:), 1, size(trig,1));
trialinds(:,trialinds(1,:)<=0) = []; trialinds(:,trialinds(end,:) > length(Accel)) = [];
Snip = Accel(floor(trialinds));
Snip = mean(Snip,2);
inds = inds/1000;
end