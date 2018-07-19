% load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataUbi.mat')
% SL = a.RemoveTwitch(SL); SL = u.AppendThreshold(SL); SL = a.AppendRT(SL); SL = a.AppendNormalizedDelay(SL); UbiSL = SL; 
% load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataIgor.mat')
% SL = a.RemoveTwitch(SL); SL = u.AppendThreshold(SL); SL = a.AppendRT(SL); SL = a.AppendNormalizedDelay(SL); IgorSL = SL;
% load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataKato.mat')
% SL = a.RemoveTwitch(SL); SL = u.AppendThreshold(SL); SL = a.AppendRT(SL); SL = a.AppendNormalizedDelay(SL); KatoSL = SL;
% 

delay = []; rt = [];
figure;
for S = 1
    switch S
        case 1
            SL = UbiSL;
        case 2
            SL = IgorSL;
        case 3
            SL = KatoSL;
    end
    for i = 1:length(SL)
        if(length(SL(i).Condition)<6 || ~strcmp(SL(i).Condition(1:6),'Contra') || isempty(SL(i).trig1) ||~isempty(SL(i).Bad) ...
                || SL(i).Condition(end)=='R')
            continue;
        end

        D = str2num(SL(i).Date); 
        Bad = [20120515,20170518];%[20120508,20120515,20120612,20120814,20120512,20111221,20120511,20120512,20120607,20120914];
        if( any(Bad==D))
            continue;
        end
        
            
        ind_var = SL(i).NormDelay;
        
        
        if(strcmp(SL(i).StimHemi,'L'))
            contRT = SL(i).rts_r;
            ipsiRT = SL(i).rts_l;
            contT = SL(i).righttrials;
            ipsiT = SL(i).lefttrials;
        elseif(strcmp(SL(i).StimHemi,'R'))
            contRT = SL(i).rts_l;
            ipsiRT = SL(i).rts_r;
            contT = SL(i).lefttrials;
            ipsiT = SL(i).righttrials;
        else
            continue;
        end
        
        disp(SL(i).Date);
        
        preI = contT(:,2)<SL(i).trig1(1);
        postI = contT(:,1)>SL(i).trig1(end);
        condI = ~(postI|preI);
        
        delay(end+1) = ind_var;
        
        pre = nanmedian(contRT(preI));
        cond = nanmedian(contRT(condI));
        post = nanmedian(contRT(postI));
        
        D = char(SL(i).Date);
        subplot(3,2,1); hold on;
        scatter(ind_var,cond-pre,'k');% text(ind_var,cond-pre,D,'FontSize',7);
        
        subplot(3,2,3); hold on;
        scatter(ind_var,post-pre,'k');% text(ind_var,post-pre,D,'FontSize',7);
        
        rt(end+1,1) = cond-pre;
        rt(end,3) = post-pre;
        
        preI = ipsiT(:,2)<SL(i).trig1(1);
        postI = ipsiT(:,1)>SL(i).trig1(end);
        condI = ~(postI|preI);
        
        pre = nanmedian(ipsiRT(preI));
        cond = nanmedian(ipsiRT(condI));
        post = nanmedian(ipsiRT(postI));
        
        subplot(3,2,2); hold on;
        scatter(ind_var,cond-pre,'k'); %text(ind_var,cond-pre,D,'FontSize',7);
        subplot(3,2,4); hold on;
        scatter(ind_var,post-pre,'k');% text(ind_var,post-pre,D,'FontSize',7);
        
        rt(end,2) = cond-pre;
        rt(end,4) = post-pre;
        
        subplot(3,2,5); hold on;
        scatter(rt(end,1),rt(end,2),'k'); % text(rt(end,1),rt(end,2),D,'FontSize',7);
    
    end
    
    subplot(3,2,1); title('Cond - Pre (contra)'); xlabel('Norm Delay (ms)'); ylabel('RT (ms)')
    subplot(3,2,2); title('Cond - Pre (ipsi)'); xlabel('Norm Delay (ms)'); ylabel('RT (ms)')
    subplot(3,2,3); title('Post - Pre (contra)'); xlabel('Norm Delay (ms)'); ylabel('RT (ms)')
    subplot(3,2,4); title('Post - Pre (ipsi)'); xlabel('Norm Delay (ms)'); ylabel('RT (ms)')
    subplot(3,2,5); title('Contra vs Ipsi (Cond-Pre)'); xlabel('Contra RT (ms)'); ylabel('Ipsi RT (ms)')
%         legend('Ubi k','Igor b','Kato r')
end

tbl = table(delay',rt(:,1),rt(:,2),rt(:,3),rt(:,4),'VariableNames',{'Delay','RT1','RT2','RT3','RT4'});
fit1 = fitlm(tbl,'RT1~Delay','RobustOpts','on');
fit2 = fitlm(tbl,'RT2~Delay','RobustOpts','on');
fit3 = fitlm(tbl,'RT3~Delay','RobustOpts','on');
fit4 = fitlm(tbl,'RT4~Delay','RobustOpts','on');
fit = fitlm(tbl,'RT2~RT1','RobustOpts','on');

subplot(3,2,6); axis off;
text(0, .95, 'Slope of LinReg, P-value', 'fontsize', 11, 'fontweight', 'bold')
text(0, .56,  [num2str(fit1.Coefficients{2,1}),'     ',num2str(fit1.Coefficients{2,4})])
text(0, .48, [num2str(fit2.Coefficients{2,1}),'     ',num2str(fit2.Coefficients{2,4})])
text(0, .40, [num2str(fit3.Coefficients{2,1}),'     ',num2str(fit3.Coefficients{2,4})])
text(0, .32, [num2str(fit4.Coefficients{2,1}),'     ',num2str(fit4.Coefficients{2,4})])
text(0, .24, [num2str(fit.Coefficients{2,1}),'     ',num2str(fit.Coefficients{2,4})])

