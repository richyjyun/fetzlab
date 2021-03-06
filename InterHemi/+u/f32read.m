function fdata=f32read(filename)
% Read the file .f32 created by trainalign function wrote by Larry. 
% The file contains the behavioral data store in the .i16 file generated by the behavioral interface train.m
% and alligned with the corresponded .bin file generated by the gugers which contains the neural signals.

% The output fdata has 9 columns:
% fdata(:,1) is the first cursor channel (Ain-0). Left Z (Green tape with black strip). values are in millivolts.
% fdata(:,2) is the first cursor channel (Ain-1). Left X. (Red tape with black strip). values are in millivolts.
% fdata(:,3) is the first cursor channel (Ain-2). Left Y. (Yellow tape with black strip).values are in millivolts.
% fdata(:,4) is the first cursor channel (Ain-3). Right Z (Green tape).values are in millivolts.
% fdata(:,5) is the first cursor channel (Ain-4). Right X (Red tape).values are in millivolts.
% fdata(:,6) is the first cursor channel (Ain-5). Right Y (Yellow tape). values are in millivolts.
% fdata(:,7) is the Trigger alignment channel which should match the guger trigger channels.
% fdata(:,8) is unused.
% fdata(:,9) is the digital Train behavioral signal channel.

% Matlab code for reading the channel 9: 
% Extract Cursor in target flag:  > inflag = bitand(fdata(:, 9), 1);
%  
% Extract target on screen flag:  > onscreen = bitand(fdata(:, 9), 2);
%  
% Extract in trial flag:  > intrial = bitand(fdata(:, 9), 4);
% 
% Extract feeder flag:  > feeder = bitand(fdata(:, 9), 8);
%  
% Extract correct response flag:  > onscreen = bitand(fdata(:, 9), 32);
%  
% Extract Box ID code:  > boxid = floor(fdata(:, 9) / 1024);
%  
% Once you have a flag extracted, you can convert it into timestamps:
% > inflag_timestamp = find((inflag(1:end-1) == 0) & (inflag(2:end) ~= 0)) + 1;

fid = fopen(filename, 'r');
fdata = fread(fid, [9, inf], 'single')'; % note the transpose to put data in columns
fclose(fid);