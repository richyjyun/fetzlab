function SL = FixTriggerDelay(SL)


for i = 1:length(SL)
    if(isempty(SL(i).trig1) || ~isempty(SL(i).Bad) || strcmp(SL(i).Condition,'nostim') || ...
            strcmp(SL(i).Condition,'tonic') || strcmp(SL(i).Condition,'Control') || strcmp(SL(i).Condition,'NaN'))
        continue;
    end
    D = SL(i).Date;
    S = SL(i).Session_Guger_Train;
    
    if(str2num(D) < 20170406)
        continue;
    end
    
    disp(SL(i).Date);
    
    
    %% Code from u.trainalign3
    
    fileprefix = [char(D),'_',char(S(2))];
    
    infile = [fileprefix '.f32'];
    datfile = [fileprefix '.bin'];
    cfgfile = [fileprefix '.cfg'];
    
    datfid = fopen(datfile, 'r');
    if (datfid < 0)
        fclose(infid);
        disp(['Error: ' fnname ' could not find daq data file ' datfile]);
        error_code = 2;
        return
    end
    infid = fopen(infile, 'r');
    if (infid < 0)
        disp(['Error: ' fnname ' could not find training data file' infile]);
        error_code = 1;
        return
    end
    
    trig_index = 1; % trying to check trig1
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Read config file
    disp('Reading guger trigs...')
    load(cfgfile,'-mat');
    nchans = length(find(UI.ch_enabled))+length(find(UI.ga_trigger)); % number of channels
    if isfield(UI, 'wd_channel')
        nchans = nchans + 1; % Window discriminator indicates daqdiscrim file.
    else
        nchans = nchans + 2; % Assume daqbinmanual was used.
    end
    chnm = cell(nchans,1); % name of channels
    chgu = zeros(nchans,1); % gUSBamp that channel was recorded on
    ii = 0;
    for iga = 1:length(UI.ga_trigger)
        ind = find(UI.ch_enabled(iga,:));
        chnm(ii+1:ii+length(ind)) = UI.ch_name(iga,ind)';
        chgu(ii+1:ii+length(ind)) = iga*ones(length(ind),1);
        ii = ii + length(ind);
        if UI.ga_trigger(iga), chnm{ii+1} = ['trig ' num2str(iga)]; chgu(ii+1) = iga; ii = ii + 1; end
    end
    if isfield(UI, 'wd_channel')
        chnm{end} = 'Discrim';  % One extra channel in daqdiscrim files.
    else
        chnm{end-1} = 'Behave1'; % Two extra chanels in daqbimanual files.
        chnm{end} = 'Behave2';
    end
    itrigch = strmatch('trig',chnm);
    fs = UI.samprate;
    
    if length(itrigch) < trig_index
        fclose(infid);
        fclose(datfid);
        disp(['Error: ' fnname ' could not find specified trigger channel in ' cfgfile]);
        error_code = 5;
        return
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Read triggers from .bin file
    % datfid = fopen(datfile, 'r');
    
    trig_offset = (itrigch(trig_index)-1) * 4;
    fseek(datfid, trig_offset, -1);
    skip_bytes = (nchans - 1) * 4;
    [samples total_samples] = fread(datfid, '*single', skip_bytes);
    fclose(datfid);
    trig1 = find(samples > 0);       % Find triggers
    intervals1 = diff([0; trig1]);   % Calcualte intervals
    trig1 = (trig1(intervals1 > 10) - 1) * 1000 / fs; % Debounce trigger and remove doublets. Time base of zero.
    intervals1 = diff([0; trig1]); % Intervals in milliseconds.
    trig1_mean = mean(intervals1(2:end));
    disp(['Mean data trigger interval = ' num2str(trig1_mean) ' ms (N=' num2str(length(trig1)) ')']);
    if (length(trig1) < 10)
        disp(['Error: found less than 10 triggers in ' datfile]);
        fclose(infid);
        error_code = 7;
        return;
    end
    
    SL(i).trig1 = trig1;
    
end

