function MakePacket(SL,fname,controlfig, split)

if ~exist('controlfig', 'var'), controlfig = figure; end
if ~exist('split', 'var'), split = 0; end

% Show first plot of behavior over all time
a.ShowMovementTimesStats(SL); 
print(gcf, '-dpsc2', fname, '-append')
close gcf

% Plot reaction times vs stimulation delay
[c,p,g,l] = u.PlotRT(SL, split, controlfig);
print(c, '-dpsc2', fname, '-append')
close(c)
print(p, '-dpsc2', fname, '-append')
close(p)
print(g, '-dpsc2', fname, '-append')
close(g)
print(l, '-dpsc2', fname, '-append')
close(l)
%print(c2, '-dpsc2', fname, '-append')
%close(c2)

for i = 1:length(SL)
    
    if(strcmp(SL(i).Bad,'1') || isempty(SL(i).trig1) || strcmp(SL(i).Condition,'nostim') || strcmp(SL(i).Condition,'tonic') || strcmp(SL(i).Condition(end),'R'))
        continue;
    end
    
    disp(char(SL(i).Date));
    
    % Plot train experiment data
    u.PlotExperiment(SL(i));
    print(gcf, '-dpsc2', fname, '-append')
    close gcf
    
    if(~strcmp(SL(i).StimHemi,'NaN'))
        % Plot neural data
        [left, right, trig, ccep] = u.PlotNeuralERP(SL(i));
        print(left, '-dpsc2', fname, '-append')
        close(left)
        print(right, '-dpsc2', fname, '-append')
        close(right)
        print(trig, '-dpsc2', fname, '-append')
        close(trig)
        %     if(strcmp(SL(i).Animal,'Ubi'))
        print(ccep, '-dpsc2', fname, '-append')
        %     end
        close(ccep)
%         SL(i) = newSL;
    end
    
end

% if(strcmp(SL(1).Animal,'Ubi'))
%     save('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataUbi.mat','SL','-v7.3');
%     [C,I] = u.PlotERPDiff(SL);
%     print(C, '-dpsc2', fname, '-append')
%     close(C)
%     print(I, '-dpsc2', fname, '-append')
%     close(I)
% end

end