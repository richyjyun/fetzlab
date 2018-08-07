function result = goodChan(blockname)
% takes in a blockname from the user and returns a numeric array of good 
% channels from the first instance of that blockname in the google
% spreadsheet.
    spdsht = GetGoogleSpreadsheet('1WLfx_3Zq1MdA2T0S6-LUTU0QqMe68vDLM8caA5_lEMc');
    token = allwords(blockname,'-');
    token = ['20',token{2}];
    for i = 1:length(spdsht)
        if (strcmp(token,spdsht{i,2}))
            GOOD_CHANNELS = spdsht{i,5};
            break
        end 
    end
    GOOD_CHANNELS = allwords(GOOD_CHANNELS,'/');
    for i = 1:length(GOOD_CHANNELS)
        result(i) = str2num(GOOD_CHANNELS{i});
end

