
ContraRT = []; IpsiRT = []; colors = []; symbols = [];

for S = 1:3
    switch S
        case 1
            SL = UbiSL;
            cs = [.5,.8,.5];
            
        case 2
            SL = IgorSL;
            cs = [.8,.5,.5];
        case 3
            SL = KatoSL;
            cs = [.5,.5,.8];
            
    end
    for i = 1:length(SL)
        cond = SL(i).Condition;
        if(strcmp(cond, 'Control') || strcmp(cond, 'nostim') || ( strcmp(cond, 'NaN') && strcmp(SL(1).Animal,'Kato')) )
            achar = '*';
        elseif ( strcmp(cond, 'tonic'))
            achar = '.';
        elseif( strcmp(cond(end),'R'))
            achar = '+';
        else
            continue;
        end
        
        colors(end+1,:) = cs;
        symbols(end+1) = achar;
        
        if(strcmp(SL(i).StimHemi,'L')|| strcmp(SL(i).StimHemi,'NaN'))
            CTrials = SL(i).righttrials;
            CRT = SL(i).rts_r;
            ITrials = SL(i).lefttrials;
            IRT = SL(i).rts_l;
        else
            ITrials = SL(i).righttrials;
            IRT = SL(i).rts_r;
            CTrials = SL(i).lefttrials;
            CRT = SL(i).rts_l;
        end
        
        trig = [];
        if(strcmp(cond(end),'R'))
            trig = SL(i).trig1;
        else
            
            
            
        if strcmp(SL(i).Animal, 'Kato')
            trig(1) =10*60*SL(i).fs;
            trig(2) = 30*60*SL(i).fs;
        elseif strcmp(SL(i).Animal, 'Igor')
            trig(1) = 45*60*1000;
            trig(2) = 75*60*1000;
        else
            trig(1) = 25*60*SL(i).fs;
            trig(2) = 60*60*SL(i).fs;
        end            
            
%             Trials = sort([SL(i).righttrials;SL(i).lefttrials]);
%             if(strcmp(SL(i).Animal,'Ubi')) % 500-1000-500
%                 trig(1) = Trials(floor(length(Trials)/4),1);
%                 trig(2) = Trials(floor(3*length(Trials)/4),1);
%             elseif(strcmp(SL(i).Animal,'Igor')) % 1300-600-1300
%                 trig(1) = Trials(floor(13*length(Trials)/32),1);
%                 trig(2) = Trials(floor(19*length(Trials)/32),1);
%             elseif(strcmp(SL(i).Animal,'Kato'))
%                 trig(1) = Trials(floor(length(Trials)/3),1);
%                 trig(2) = Trials(floor(2*length(Trials)/3),1);
%             end   
        end
        
        bounds = [find(CTrials(:,2)<trig(1),1,'last'),find(CTrials(:,1)>trig(end),1)];
        
        ContraRT(end+1,1) = nanmedian(CRT(1:bounds(1)));
        ContraRT(end,2) = nanmedian(CRT(bounds(1)+1:bounds(2)-1));
        if(isnan(ContraRT(end,2)))
            keyboard;
        end
        ContraRT(end,3) = nanmedian(CRT(bounds(2):end));
        
        bounds = [find(ITrials(:,2)<trig(1),1,'last'),find(ITrials(:,1)>trig(end),1)];
        
        IpsiRT(end+1,1) = nanmedian(IRT(1:bounds(1)));
        IpsiRT(end,2) = nanmedian(IRT(bounds(1)+1:bounds(2)-1));
        IpsiRT(end,3) = nanmedian(IRT(bounds(2):end));
        
        
    end
end

figure; subplot(1,2,1); hold on;
for i = 1:size(ContraRT,1)
    hold on; scatter([1,2,3],ContraRT(i,:),75,char(symbols(i)),'MarkerEdgeColor',colors(i,:));
    hold on; plot([1,2,3],ContraRT(i,:),'Color',colors(i,:));
end
% b = bar([1,2,3],mean(ContraRT),'w');
% uistack(b,'bottom');
set(gca,'XTick',1:3); xlim([0.5,3.5]);
set(gca,'XTickLabel',{'PRE','STIM','POST'});
title('Contra Trials'); ylim([200,500]); ylabel('RT (ms)')

subplot(1,2,2); hold on;
for i = 1:size(IpsiRT,1)
    hold on; scatter([1,2,3],IpsiRT(i,:),75,char(symbols(i)),'MarkerEdgeColor',colors(i,:));
    hold on; plot([1,2,3],IpsiRT(i,:),'Color',colors(i,:));
end
% b = bar([1,2,3],mean(IpsiRT),'w');
% uistack(b,'bottom');
set(gca,'XTick',1:3); xlim([0.5,3.5]);
set(gca,'XTickLabel',{'PRE','STIM','POST'});
title('Ipsi Trials'); ylim([200,500]);

subplot(1,2,1); hold on; 
h = [];
h(1) = scatter(nan,nan,25,[.8,.5,.5],'filled');
h(2) = scatter(nan,nan,25,[.5,.5,.8],'filled');
h(3) = scatter(nan,nan,25,[.5,.8,.5],'filled');
legend(h,'Monkey I','Monkey K','Monkey U'); legend('boxoff')

subplot(1,2,2); hold on;
h = [];
h(1) = scatter(nan,nan,'k*');
h(2) = scatter(nan,nan,'k.');
h(3) = scatter(nan,nan,'k+');
legend(h,'No Stim','Tonic','Random'); legend('boxoff')

% 
% [h,p] = ttest(ContraRT(:,1),ContraRT(:,3))
% [h,p] = ttest(ContraRT(:,2),ContraRT(:,3))
% [h,p] = ttest(ContraRT(:,2),ContraRT(:,1))
% [h,p] = ttest(IpsiRT(:,2),IpsiRT(:,3))
% [h,p] = ttest(IpsiRT(:,1),IpsiRT(:,3))
% [h,p] = ttest(IpsiRT(:,2),IpsiRT(:,1))
