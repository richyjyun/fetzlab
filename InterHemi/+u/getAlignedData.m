function SL = getAlignedData(SL)


for i = 1:length(SL)
    if(isempty(SL(i).trig1) || ~isempty(SL(i).Bad) || strcmp(SL(i).Condition,'nostim') || ...
            strcmp(SL(i).Condition,'tonic') || strcmp(SL(i).Condition,'Control') || strcmp(SL(i).Condition,'NaN'))
        continue;
    end
    
    disp(SL(i).Date);
    
    % Loading in the data.
    D = SL(i).Date;
    S = SL(i).Session_Guger_Train;
    Session = [char(D),'_',char(S(2))];
    trainfile = [Session,'.f32'];
    gugfile = Session;
    if(~exist(trainfile) || ~exist([gugfile,'.bin']) || ~exist([gugfile,'.cfg']))
        continue;
    end
    
    StimChns = strsplit(SL(i).Stim_Loc,'/');
    n = char(StimChns(1));
    StimChns(2) = cellstr([n(1:4),char(StimChns(2))]);
    
    % down sampling rate
    dwn = 5;
    [accel, trig1, trig2, lefttrials, righttrials, fs, lefttrialsuccess, righttrialsuccess] = u.LoadTrainUbi(trainfile,dwn);
    
    x = 0:1/fs:1/fs*(length(accel)-1);
    newx = linspace(0,x(end),(length(accel)-1)*1000/fs);
    SL(i).accel_raw_l = interp1(x,accel(:,1),newx);
    SL(i).accel_raw_r = interp1(x,accel(:,2),newx);
    % ones after already has gotten their triggers fixed using neural data, 
    % not aligned data, due to randomizing code delay bug
    if(str2num(D) < 20170406) 
        SL(i).trig1 = round(trig1*1000/fs);
    end
    if(str2num(D) == 20170127) 
        SL(i).trig1 = SL(i).trig1(12:end);
    end
    if(str2num(D) == 20170411)
        SL(i).trig1(end) = [];
    end
    SL(i).trig2 = round(trig2*1000/fs);
    SL(i).lefttrials = round(lefttrials*1000/fs);
    SL(i).righttrials = round(righttrials*1000/fs);
    SL(i).lefttrialsuccess = lefttrialsuccess;
    SL(i).righttrialsuccess = righttrialsuccess;
    
end


end