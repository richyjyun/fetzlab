% For fixing train data (including triggers) to the aligned version 

function SL = getAlignedDataKato(SL)


for i = 1:length(SL)
    if(isempty(SL(i).trig1) || ~isempty(SL(i).Bad) || strcmp(SL(i).Condition,'nostim') || ...
            strcmp(SL(i).Condition,'tonic') || strcmp(SL(i).Condition,'Control') || strcmp(SL(i).Condition,'NaN'))
        continue;
    end
    
    disp(SL(i).Date);
    
    %Loading in data
    if(strcmp(SL(i).Date,'20150108')),continue;end;
    D = SL(i).Date;
    Sessions = strsplit(char(SL(i).Sessions),'_');
    dwn = 5;
    trig1 = []; lefttrials = []; righttrials = []; lefttrialsuccess = []; righttrialsuccess = []; data = []; last = 0; accel = [];
    trig2 = [];
    for j = 1:3
        trainfile = [char(D),'_',char(Sessions(j)),'.f32'];
        gugfile = [char(D),'_',char(Sessions(j))];
        if(~exist([gugfile,'.i16']) || ~exist([gugfile,'.bin']) || ~exist([gugfile,'.cfg']))
            return;
        end
        if(~exist(trainfile))
            disp(['making f32 file for ' gugfile]);
            cd('F:\Dropbox\repos\abogaard\efetz\DualAccelRt\code');
            err = utils.trainalign2(gugfile, 1);
        end
        
        [acc, t1, lt, rt, fs1, lts, rts] = u.LoadTrainKato(trainfile,dwn);
        chn = 2; % trig1 channel 
        % Code borrowed from trainalign
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Read config file
        load([gugfile,'.cfg'],'-mat');
        fs = UI.samprate;
        nchans = length(find(UI.ch_enabled))+length(find(UI.ga_trigger)); % number of channels
        if isfield(UI, 'wd_channel')
            nchans = nchans + 1; % Window discriminator indicates daqdiscrim file.
        else
            nchans = nchans + 2; % Assume daqbinmanual was used.
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Read bin file
        datfid = fopen([gugfile,'.bin'], 'r');
        if (datfid < 0)
            fclose(datfid);
            disp(['Error: ' fnname ' could not find daq data file ', gugfile,'.bin']);
            error_code = 2;
            return
        end
        trig_offset = (chn-1) * 4;
        fseek(datfid, trig_offset, -1);
        skip_bytes = (nchans - 1) * 4;
        [samples total_samples] = fread(datfid, '*single', skip_bytes);
        fclose(datfid);
        trig = find(samples > 0);       % Find triggers
        intervals1 = diff([0; trig]);   % Calcualte intervals
        trig = (trig(intervals1 > 10) - 1); % Debounce trigger and remove doublets. Time base of zero.
        
        lefttrials = [lefttrials;last+lt];
        righttrials = [righttrials;last+rt];
        lefttrialsuccess = [lefttrialsuccess;lts];
        righttrialsuccess = [righttrialsuccess;rts];
        trig2 = [trig2;last+trig];
        if(j == 2)
            trig1 = last+trig;
        end
        accel = [accel;acc];
        last = last+length(acc);
    end
     
    
    x = 0:1/fs1:1/fs1*(length(accel)-1);
    newx = linspace(0,x(end),(length(accel)-1)*1000/fs1);
    SL(i).accel_raw_l = interp1(x,accel(:,1),newx);
    SL(i).accel_raw_r = interp1(x,accel(:,2),newx);
    SL(i).trig1 = round(trig1*1000/fs);
    SL(i).trig2 = round(trig2*1000/fs);
    SL(i).lefttrials = round(lefttrials*1000/fs1);
    SL(i).righttrials = round(righttrials*1000/fs1);
    SL(i).lefttrialsuccess = lefttrialsuccess;
    SL(i).righttrialsuccess = righttrialsuccess;
    
end


end