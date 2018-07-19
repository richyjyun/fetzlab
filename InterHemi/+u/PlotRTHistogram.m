% figure;
left = [];
right = [];
for s = 1
%     if(s==1) 
%         SL = UbiSL;
%         M = 'Ubi';
%     elseif(s==2)
%         SL = IgorSL;
%         M = 'Igor';
%     elseif(s==3)
%         SL = KatoSL;
%         M = 'Kato';
%     end
    
    LRT = []; RRT = [];
    
    for i = 1:length(SL)
        if (~isempty(SL(i).Bad) || strcmp(SL(i).Condition, 'Control')...
                || strcmp(SL(i).Condition, 'NaN') || strcmp(SL(i).Condition, 'nostim')...
                || isempty(SL(i).trig1) || strcmp(SL(i).Condition(end),'R'))
            continue;
        end
        
        lind = find(SL(i).lefttrials(:,2) < SL(i).trig1(1),1,'last');
        rind = find(SL(i).righttrials(:,2) < SL(i).trig1(1),1,'last');
        left(end+1) = nanmedian(SL(i).rts_l(1:lind));
        right(end+1) = nanmedian(SL(i).rts_r(1:rind));
        LRT = [LRT;SL(i).rts_l(1:lind)];
        RRT = [RRT;SL(i).rts_r(1:rind)]; 
    end
    
%     bins = 0:10:600;
%     h1 = subplot(3,2,(s-1)*2+1); histogram(LRT,bins,'FaceColor','k','Edgecolor','none','FaceAlpha',1); 
%     title([M,' Left Hand RT']);
%     h2 = subplot(3,2,(s-1)*2+2); histogram(RRT,bins,'FaceColor','k','Edgecolor','none','FaceAlpha',1); 
%     title([M,' Right Hand RT']);
%     linkaxes([h1,h2],'y');
end

%%
% figure;
left = [];
right = [];
for s = 1
%     if(s==1) 
%         SL = UbiSL;
%         M = 'Ubi';
%     elseif(s==2)
%         SL = IgorSL;
%         M = 'Igor';
%     elseif(s==3)
%         SL = KatoSL;
%         M = 'Kato';
%     end
    
    LRT = []; RRT = [];
    
    for i = 1:length(SL)
        if (~isempty(SL(i).Bad) || strcmp(SL(i).Condition, 'Control')...
                || strcmp(SL(i).Condition, 'NaN') || strcmp(SL(i).Condition, 'nostim')...
                || isempty(SL(i).trig1) || strcmp(SL(i).Condition(end),'R'))
            continue;
        end
        
        lind = find(SL(i).lefttrialsuccess == 0);
        rind = find(SL(i).righttrialsuccess == 0);
        left = [left,diff(SL(i).lefttrials(lind,:)')];
        right = [right,diff(SL(i).righttrials(rind,:)')];

    end
    
end



