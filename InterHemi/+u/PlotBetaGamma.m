% function PlotBetaGamma(SLNeuro,SL,fname)

fname = fullfile('F:\Dropbox\repos\abogaard\efetz\U\packets', 'BetaGammaUbi.ps');
if(exist(fname))
    delete(fname);
end
normalized = 0;
epoch = 'cond';
dates = extractfield(SL,'Date'); window = 0.5;
Beta = struct('L',[],'R',[]); Gamma = Beta; RT = Beta;
for i = 1:length(SLNeuro)
    inds = SLNeuro(i).tneuro > -window & SLNeuro(i).tneuro < window;
    
    %% find proper SL, get change in RT
    % Doing post - pre for igor, can change to cond-pre for ubi
    SLInd = find(strcmp(dates,SLNeuro(i).Date)); trig = SL(SLInd).trig1;
    D = str2num(SL(SLInd).Date);
    if(~strcmp(SL(SLInd).Condition(1:6),'Contra'))
        continue;
    end
    
    if(D >= 20170226 && D <= 20170306)
        continue;
    end
    bounds = [find(SL(SLInd).lefttrials(:,2)<trig(1),1,'last'),find(SL(SLInd).lefttrials(:,2)>trig(end),1)];
    
    if(strcmp(SL(SLInd).Animal,'Igor'))
        RT.L(end+1) = nanmedian(SL(SLInd).rts_l(bounds(2):end)) - nanmedian(SL(SLInd).rts_l(1:bounds(1)));
        bounds = [find(SL(SLInd).righttrials(:,2)<trig(1),1,'last'),find(SL(SLInd).righttrials(:,2)>trig(end),1)];
        RT.R(end+1) = nanmedian(SL(SLInd).rts_r(bounds(2):end)) - nanmedian(SL(SLInd).rts_r(1:bounds(1)));
    elseif(strcmp(SL(SLInd).Animal,'Ubi'))
        RT.L(end+1) = nanmedian(SL(SLInd).rts_l(bounds(1)+1:bounds(2)-1)) - nanmedian(SL(SLInd).rts_l(1:bounds(1)));
        bounds = [find(SL(SLInd).righttrials(:,2)<trig(1),1,'last'),find(SL(SLInd).righttrials(:,2)>trig(end),1)];
        RT.R(end+1) = nanmedian(SL(SLInd).rts_r(bounds(1)+1:bounds(2)-1)) - nanmedian(SL(SLInd).rts_r(1:bounds(1)));
    end
    
    %% using end - 1 to get difference from post for now. Can't really do
    % cond with igor, can change to 2-1 for Ubi
    if(strcmp(epoch,'post'))
        if(normalized)
            base = round(SLNeuro(i).fs*0.1);
            a = SLNeuro(i).Lbeta(end,:,inds); b = SLNeuro(i).Lbeta(1,:,inds);
            Beta.L(end+1,:) = (min(a,[],3)-min(b,[],3))./(mean(a(1,:,1:base),3)-mean(b(1,:,1:base),3));
            a = SLNeuro(i).Rbeta(end,:,inds); b = SLNeuro(i).Rbeta(1,:,inds);
            Beta.R(end+1,:) = (min(a,[],3)-min(b,[],3))./(mean(a(1,:,1:base),3)-mean(b(1,:,1:base),3));
            a = SLNeuro(i).Lgamma(end,:,inds); b = SLNeuro(i).Lgamma(1,:,inds);
            Gamma.L(end+1,:) = (max(a,[],3)-max(b,[],3))./(mean(a(1,:,1:base),3)-mean(b(1,:,1:base),3));
            a = SLNeuro(i).Rgamma(end,:,inds); b = SLNeuro(i).Rgamma(1,:,inds);
            Gamma.R(end+1,:) = (max(a,[],3)-max(b,[],3))./(mean(a(1,:,1:base),3)-mean(b(1,:,1:base),3));
        else
            Beta.L(end+1,:) = min(SLNeuro(i).Lbeta(end,:,inds),[],3)-min(SLNeuro(i).Lbeta(1,:,inds),[],3);
            Beta.R(end+1,:) = min(SLNeuro(i).Rbeta(end,:,inds),[],3)-min(SLNeuro(i).Rbeta(1,:,inds),[],3);
            Gamma.L(end+1,:) = max(SLNeuro(i).Lgamma(end,:,inds),[],3)-max(SLNeuro(i).Lgamma(1,:,inds),[],3);
            Gamma.R(end+1,:) = max(SLNeuro(i).Rgamma(end,:,inds),[],3)-max(SLNeuro(i).Rgamma(1,:,inds),[],3);
        end
    elseif(strcmp(epoch,'cond'))
        Beta.L(end+1,:) = min(SLNeuro(i).Lbeta(2,:,inds),[],3)-min(SLNeuro(i).Lbeta(1,:,inds),[],3);
        Beta.R(end+1,:) = min(SLNeuro(i).Rbeta(2,:,inds),[],3)-min(SLNeuro(i).Rbeta(1,:,inds),[],3);
        Gamma.L(end+1,:) = max(SLNeuro(i).Lgamma(2,:,inds),[],3)-max(SLNeuro(i).Lgamma(1,:,inds),[],3);
        Gamma.R(end+1,:) = max(SLNeuro(i).Rgamma(2,:,inds),[],3)-max(SLNeuro(i).Rgamma(1,:,inds),[],3);
    end
    
    if(strcmp(SL(SLInd).Animal,'Ubi'))
        if(str2num(SL(SLInd).Stim_Loc(5)) == 3) % remove those channels
            Beta.L(end,[21,22]) = nan; Beta.R(end,[21,22]) = nan;
            Gamma.L(end,[21,22]) = nan; Gamma.R(end,[21,22]) = nan;
        elseif(str2num(SL(SLInd).Stim_Loc(5)) == 4)
            Beta.L(end,[23,24]) = nan; Beta.R(end,[23,24]) = nan;
            Gamma.L(end,[23,24]) = nan; Gamma.R(end,[23,24]) = nan;
        end
    end
        
    
end

% plotAll(Beta,Gamma,RT,1:16)
% scatterAll(Beta,Gamma,RT,1:9)
% scatterAll(Beta,Gamma,RT,10:18)


f = plotfreqs(SLNeuro,Beta,RT);
for t = 1:length(f)
    set(0, 'CurrentFigure', f(t))
    a = axes; t1 = title('\Delta Beta vs \Delta RT');
    a.Visible = 'off'; t1.Visible = 'on';
    print(f(t), '-dpsc2', fname, '-append')
    %     close(f(t))
end
f = plotfreqs(SLNeuro,Gamma,RT);
for t = 1:length(f)
    set(0, 'CurrentFigure', f(t))
    a = axes; t1 = title('\Delta Gamma vs \Delta RT');
    a.Visible = 'off'; t1.Visible = 'on';
    print(f(t), '-dpsc2', fname, '-append')
    %     close(f(t))
end

function f = plotfreqs(SLNeuro,freq,RT)

f = []; f(1) = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]); %set(f(1),'visible','off');
for chn = 1:length(SLNeuro(1).chID)
    fInd = chn*2-1;
    if(fInd > 36)
        if(length(f)<2)
            f(2) = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]); %set(f(1),'visible','off');
        end
        fInd = chn*2-1-36;
    end
    
    nbins = 8;
    
    subaxis(6, 6, fInd, 'padding', .01, 'spacing', .01, 'spacingvert', .03);
    scatter(RT.L,freq.L(:,chn),'k.');% ylim([-1,1]);  %ylim([-0.3,0.3])
%     [~,edges] = histcounts(RT.L,nbins); Y = discretize(RT.L,edges);
%     binned = zeros(length(edges)-1,1); 
%     for bin = 1:nbins
%         binned(bin) = nanmean(freq.L(Y==bin,chn));
%     end
%     mid = edges(1:end-1) + diff(edges)/2;
%     bar(mid,binned,1); xticks(sort([edges(1),0,edges(end)])); xlim([edges(1),edges(end)]);    
%     title([SLNeuro(1).chnm(SLNeuro(1).chID(chn)),' Left'],'fontsize',7);
    
    %     tbl = table(RT.L',freq.L(:,chn),'VariableNames',{'RT','Freq'});
    %     linfit = fitlm(tbl,'Freq~RT','RobustOpts','on');
    %     h = plot(linfit); h(1).Marker = '.'; h(1).Color = [0,0,0]; legend HIDE;
    %     title([SLNeuro(chn).chnm(SLNeuro(chn).chID(chn)),' Left;',...
    %         ' Slope: ',num2str(linfit.Coefficients{2,1}),'; P-value: ',num2str(linfit.Coefficients{2,2})],'fontsize',7);
    %     xlabel(''); ylabel('');   set(gca,'fontsize',7);
    
    subaxis(6, 6, fInd+1, 'padding', .01, 'spacing', .01, 'spacingvert', .03);
    scatter(RT.R,freq.R(:,chn),'r.');%% ylim([-1,1]);  %ylim([-0.3,0.3])
%     [~,edges] = histcounts(RT.R,nbins); Y = discretize(RT.R,edges);
%     binned = zeros(length(edges)-1,1); 
%     for bin = 1:nbins
%         binned(bin) = nanmean(freq.R(Y==bin,chn));
%     end
%     mid = edges(1:end-1) + diff(edges)/2;
%     bar(mid,binned,1); xticks(sort([edges(1),0,edges(end)])); xlim([edges(1),edges(end)]);
%     title([SLNeuro(chn).chnm(SLNeuro(chn).chID(chn)),' Right'],'fontsize',7);
    
    %     tbl = table(RT.L',freq.R(:,chn),'VariableNames',{'RT','Freq'});
    %     linfit = fitlm(tbl,'Freq~RT','RobustOpts','on');
    %     h = plot(linfit); h(1).Marker = '.'; h(1).Color = [1,0,0]; legend HIDE;
    %     title([SLNeuro(chn).chnm(SLNeuro(chn).chID(chn)),' Right;',...
    %         ' Slope: ',num2str(linfit.Coefficients{2,1}),'; P-value: ',num2str(linfit.Coefficients{2,2})],'fontsize',7);
    %     xlabel(''); ylabel(''); set(gca,'fontsize',7);
    
end

end

function f = plotAll(Beta,Gamma,RT,chns)
nbins = 10;
f = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]); %set(f(1),'visible','off');
[~,Ledges] = histcounts(RT.L,nbins); Lbins = discretize(RT.L,Ledges);
[~,Redges] = histcounts(RT.R,nbins); Rbins = discretize(RT.R,Redges);
BetaAll = cell(2,6);  GammaAll = cell(2,6);

for bin = 1:nbins
    BetaAll{1,bin} = reshape(Beta.L(Lbins==bin,chns),1,size(Beta.L(Lbins==bin,chns),1)*size(Beta.L(Lbins==bin,chns),2));
    BetaAll{2,bin} = reshape(Beta.R(Rbins==bin,chns),1,size(Beta.R(Rbins==bin,chns),1)*size(Beta.R(Rbins==bin,chns),2));
    GammaAll{1,bin} = reshape(Gamma.L(Lbins==bin,chns),1,size(Gamma.L(Lbins==bin,chns),1)*size(Gamma.L(Lbins==bin,chns),2));
    GammaAll{2,bin} = reshape(Gamma.R(Rbins==bin,chns),1,size(Gamma.R(Rbins==bin,chns),1)*size(Gamma.R(Rbins==bin,chns),2));
end

Lmid = Ledges(1:end-1) + diff(Ledges)/2;
subplot(2,2,1); bar(Lmid,cellfun(@nanmean,BetaAll(1,:)),1); title('Beta Left');
subplot(2,2,3); bar(Lmid,cellfun(@nanmean,GammaAll(1,:)),1); title('Gamma Left');

Rmid = Redges(1:end-1) + diff(Redges)/2;
subplot(2,2,2); bar(Rmid,cellfun(@nanmean,BetaAll(2,:)),1); title('Beta Right');
subplot(2,2,4); bar(Rmid,cellfun(@nanmean,GammaAll(2,:)),1); title('Gamma Right');
end


function f = scatterAll(Beta,Gamma,RT,chns)
f = figure;

for chn = 1:length(chns)
    subplot(2,2,1); hold on; scatter(RT.L,Beta.L(:,chns(chn)),'k.');
    subplot(2,2,2); hold on; scatter(RT.R,Beta.R(:,chns(chn)),'k.');
    subplot(2,2,3); hold on; scatter(RT.L,Gamma.L(:,chns(chn)),'k.');
    subplot(2,2,4); hold on; scatter(RT.R,Gamma.R(:,chns(chn)),'k.');
end

subplot(2,2,1); title('Beta Left');   ylim([-1,1])
subplot(2,2,2);  title('Beta Right'); ylim([-1,1])
subplot(2,2,3);  title('Gamma Left'); ylim([-0.5,0.5])
subplot(2,2,4); title('Gamma Right'); ylim([-0.5,0.5])

% BetaLAll = Beta.L(:,chns); BetaLAll = reshape(BetaLAll',numel(BetaLAll),1);
% GammaLAll = Beta.L(:,chns); GammaLAll = reshape(GammaLAll',numel(GammaLAll),1);
% 
% RTL = repmat(RT.L',length(chns),1);
% tbl = table(RTL,GammaLAll,'VariableNames',{'RT','Gamma'});
% linfit = fitlm(tbl,'Gamma~RT','RobustOpts','on');

end


