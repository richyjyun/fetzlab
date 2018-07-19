function SL = RemoveTwitch(SL,dwn)

if ~exist('dwn', 'var'), dwn=1; end

for i = 1:length(SL)
    
    % Code for blanking out twitches
    Accel = [];
    if(strcmp(SL(i).Animal,'Ubi'))
        if(strcmp(SL(i).Condition(end),'0'))
            blank = [25,125];
        elseif(strcmp(SL(i).Condition(end),'2'))
            blank = [25,75]; %ms after stim to blank out
        else
            continue;
        end
    elseif(strcmp(SL(i).Animal,'Igor'))
        if(strcmp(SL(i).Condition,'Contra_0'))
            blank = [0,200];
        elseif(strcmp(SL(i).Condition,'Contra_250'))
            blank = [50,125];
        else
            continue;
        end
    elseif(strcmp(SL(i).Animal,'Kato'))
        if(strcmp(SL(i).Condition(end),'2'))
            blank = [125,175]; %ms after stim to blank out
        else
            continue;
        end
    end
    if(any(strmatch('Contra',SL(i).Condition)) && ~strcmp(SL(i).Condition,'Control') && ~isempty(SL(i).trig1) && ~strcmp(SL(i).Condition(end),'R'))
        if(strcmp(SL(i).StimHemi,'L'))
            Accel = SL(i).accel_raw_r;
        elseif(strcmp(SL(i).StimHemi,'R'))
            Accel = SL(i).accel_raw_l;
        end
    end
    
    if(~isempty(Accel))
        %         if(any(strmatch('Ipsi',SL(i).Condition)))
        %             keyboard;
        %         end
        for j = 1:length(SL(i).trig1)
            ind = floor(SL(i).trig1(j)*SL(i).fs/1000) + str2num(SL(i).Stim_Delay);
            Accel(ind+blank(1):ind+blank(2)) = linspace(Accel(ind+blank(1)-1),Accel(ind+blank(2)+1),blank(2)-blank(1)+1);
        end
        if(strcmp(SL(i).StimHemi,'L'))
            SL(i).accel_raw_r = Accel;
        elseif(strcmp(SL(i).StimHemi,'R'))
            SL(i).accel_raw_l = Accel;
        end
    end
    % End of blanking out code
    
end

end