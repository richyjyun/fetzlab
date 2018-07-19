function SL = RemoveTwitch2(SL,dwn)

if ~exist('dwn', 'var'), dwn=1; end

for i = 1:length(SL)
    
    % Code for blanking out twitches
    Accel = [];
    if(strcmp(SL(i).Animal,'Ubi'))
        if(strcmp(char(SL(i).Date),'20170411'))
            blank = [25,100];
        else
            blank = [50,150]; %ms after stim to blank out
        end
    elseif(strcmp(SL(i).Animal,'Igor'))
        blank = [10,40];
    end
    if(~strcmp(SL(i).Animal,'Kato') && any(strmatch('Contra',SL(i).Condition)) && ~strcmp(SL(i).Condition,'Control') && ~isempty(SL(i).trig1) && ~strcmp(SL(i).Condition(end),'R'))
        if(strcmp(SL(i).StimHemi,'L'))
            Accel = SL(i).accel_raw_r;
        elseif(strcmp(SL(i).StimHemi,'R'))
            Accel = SL(i).accel_raw_l;
        end
    end
    
    blank = blank/1000*SL(i).fs*dwn;
    if(~isempty(Accel))
        %         if(any(strmatch('Ipsi',SL(i).Condition)))
        %             keyboard;
        %         end
        for j = 1:length(SL(i).trig1)
            ind = floor(SL(i).trig1(j) + str2num(SL(i).Stim_Delay)/1000*SL.fs*dwn);
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