function varargout = Neurochip3Gui(varargin) %#ok<FNDEF>
varargout = {};

% GUI for Neurochip3
% nc307 -- Stimulation safety limits placed on total number stimulations
%   received inside of a specified time window.
% nc320 -- Separate refractory period for window discriminators.
%   Event selection: Off, Window Discrim, Periodic, Discrim + Periodic
%   Support for 8 Conditions, and 6 stimulus channels. Changed the GUI
%   numbering of channels/conditions etc to start at index 1. Safety limits
%   are now on every stimulus channel.
% nc321 -- Rearranged Stim tab.  Intan Fast Settle has more options.
% nc322 -- Event refractory period working with sliding window.
%   Stimulaton requirements count>N, count<=N working with sliding window.
% nc323 -- Improvements in filter plots.  Intan Aux3 stim pin selection.
%   Z Check tab for impedance testing.  Impedances seem fairly accurate
%   from 20kOhm to 1MOhm.  Neurochip front end capacitance and resistance
%   gives a minumum impedance of about 8kOhm at 1kHz.

ui.version = 323;        % GUI version
ui.name = ['nc' num2str(ui.version)];
ui.path = 'C:\nc3\';     % Where to store settings files
ui.datapath = 'C:\data'; % Where to store downloaded data
ui.sddrive = 'F';        % Default drive letter for SDCard reader.
ui.last_settings = [ui.path 'NC320_Last_Settings.txt'];
ui.session_id = 0;       % Session identifier.
ui.show_10MOhm_scale = 0;% 0 don't show 10MOhm impedance test scale, 1 = show.  This scale can be very noisey and doesn't seem accurate.
ui.show_test_tab = 0;    % 0 don't show Test tab, 1 show test tab.
ui.show_test_button = 0; % 0 don't show TestButton, 1 show the test button.
ui.show_remote_channels = 0; % 0 don't show remote channels, 1 show remote channels pulldown menu (for a future addition)
ui.save_in_progress = 0; % 1 while the ncsave mex file is running.

% Conversion factors and constants

ui.ms2fs = 20; % Converts milliseconds to 20kHz A/D sample rate (fs).
ui.adu2uv = 0.195;         % From Intan data sheet.  Converts A/D units to microvolts.
ui.uv2adu = 1 / ui.adu2uv; % Converts microvolts to A/D units.
ui.accel2adu = 2^14 / 1000;      % 16834 accel units is 1g, but we give levels in increments of 0.001g
ui.accelTemp2adu = 64;         % 1 degrees celsius = 64 digital temp units.
ui.fs = 20000; % Max A/D sampling rate
ui.fs2 = ui.fs / 2;
ui.fs2sec = 1 / ui.fs;

% Storage for sending and receiving data from neurochip.

ui.registerBlockSize = 2048;  % Size of GUI->CPU parameter block in bytes.
ui.r = zeros(ui.registerBlockSize,1,'uint8');

% Parameters for Neurochip3
% These are values that are stored in parameter files.

p.subject = 'NC3';       % Prefix name for settings files and downloads.
p.version = ui.version;  % Parameter set version
p.session_id = ui.session_id;  % Session when parameters were saved
p.neurochip_firmware = 'Not_Loaded';  % Last downloaded neurochip firmware version string.
p.date = date;
p.file_name = '';        % Original file name. ^^^TODO make sure to save this
p.analog_lower_bandwidth = 10;      % 10 Hz
p.analog_lower_bandwidth_code = 12; % Intan code for low 10 Hz cut
p.analog_upper_bandwidth = 7500;    % 7500 hz
p.analog_upper_bandwidth_code = 3;  % Intan code for high 7500 Hz cut
p.dsp_high_pass_code = 0;           % Intan code for the high pass dsp off
p.dsp_high_pass_info = 'Off';
ui.channel_names = {'Ch1';'Ch2';'Ch3';'Ch4';'Ch5';'Ch6';'Ch7';'Ch8'; ...
    'Ch9';'Ch10';'Ch11';'Ch12';'Ch13';'Ch14';'Ch15';'Ch16';...
    'Ch17';'Ch18';'Ch19';'Ch20';'Ch21';'Ch22';'Ch23';'Ch24';...
    'Ch25';'Ch26';'Ch27';'Ch28';'Ch29';'Ch30';'Ch31';'Ch32'; ...
    'Aux1';'Aux2';'Aux3';'AccelX';'AccelY';'AccelZ';'AccelTemp';'AccelMagnitude'; ...
    'Rem1';'Rem2';'Rem3';'Rem4';'Rem5';'Rem6';'Rem7';'Rem8'};
ui.channel_Ch_index = 1:32;
ui.channel_Aux_index = 33:35;
ui.channel_Accel_index = 36:40;
ui.channel_Rem_index = 41:48;
ui.channel_count = length(ui.channel_names);
p.channel_names = ui.channel_names; % Default channel names.
p.channel_count = ui.channel_count;
p.intan_channels = 32;  % Number of Intan channels: 16 or 32
p.remote_channels = 0;  % Number of remote emg channels.
p.channel_rate = repmat(5000, p.channel_count, 1);  % Base sampling rate for Intan channels
p.channel_rate(ui.channel_Aux_index) = 0; % Sample rate of remote AUX channels.
p.channel_rate(ui.channel_Accel_index) = 100; % Sample rate of acceleromeer channels.
p.channel_rate(ui.channel_Rem_index) = 2000; % Sample rate of remote EMG channels.
p.event_count = 8; % Number of event generators defined
p.event_method = zeros(p.event_count, 1);    % 0 = 0ff, 1 = Window Discriminator, 2 = Periodic, 3 = Discrim + Perioic.
p.event_source1 = (1:p.event_count)';        % 0 = none, 1 = first analog channel...
p.event_source2 = zeros(p.event_count, 1);   % 0 = none, 1 = first analog channel...
p.event_transform = zeros(p.event_count, 1); % 0 = none, 1 = rectify.

p.event_lower_bandwidth = repmat(10, p.event_count, 1);   % The event filter is a 1st order butterworth bandpass.
p.event_upper_bandwidth = repmat(7500, p.event_count, 1); % It uses 24-bit coefficients sent to the FPGA.
p.event_filter_b = repmat([11847856 0 -11847856], p.event_count, 1); % Calculated filter coefficients for butter(1,[10 7500]/10000)
p.event_filter_a = repmat([16777216 -9784230 -6918495], p.event_count, 1);

p.event_ratio = ones(p.event_count, 1);         % 1..N, Periodic event decimator to allow for more normal like distributions
p.event_interval_min = repmat(5, p.event_count, 1);  % Minimum periodic interval.
p.event_interval_max = zeros(p.event_count, 1);  % Maximum periodic interval. 0 = off
p.event_interval_exp = zeros(p.event_count, 1);  % Exponential distribution target mean (Hz).  0 = will default to uniform distribution.

p.event_threshold = repmat(100, p.event_count, 1);    % Initial detection threshold in uV.
p.event_window1_delay = repmat(5, p.event_count, 1);  % Millisecond delay from threshold detection
p.event_window1_max = repmat(100, p.event_count, 1);  % Signal must pass >= min but < Max
p.event_window1_min = zeros(p.event_count, 1);
p.event_window2_delay = repmat(10, p.event_count, 1); % Second window
p.event_window2_max = zeros(p.event_count, 1);
p.event_window2_min = repmat(-100, p.event_count, 1);
p.event_refractory_period = repmat(5, p.event_count, 1);   % Refractory period when refractory_count reached.
p.event_refractory_count = ones(p.event_count, 1);    % Count triggering refractory period.
p.event_refractory_window = repmat(5, p.event_count, 1);   % Counts must all happen inside of given millisecond window.

p.stim_power = 0; % 0 = off, Voltage compliance = 10 * p.stim_power
p.stim_artifact_time = 0; % Milliseconds to activate artifact suppression after each stimulus.
p.stim_aux3_input = 0;    % Input selection for Intan Aux3.  0=Tissue Ground (PCal_VC0), 1=A1_VC1, ..., 6=C0_VC7, 7=Tissue Ground (NCal_VC7), >=8 Not Connected
p.stim_fast_settle = 0;   % 0..7  Number of 5kHz Intan sample packets to set the fast settle bit after each stimulation pulse starts.
p.stim_fast_settle_menu = {'Off'};
for item = 1:12
    p.stim_fast_settle_menu{item+1} = [num2str(200 * item) ' uS on each pulse'];
    % Uncomment the following line to enable post pulse options.
    %p.stim_fast_settle_menu{item+11} = [num2str(200 * item) ' uS post pulse'];
end
p.stim_cond_count = 8;    % Number of possible stimulus conditions
p.stim_cond_start = 0;    % Starting condition (0..p.stim_cond_count).  0 for Condition1, 8 for no stimulation.
p.stim_channel_count = 6; % Number of stimulus channels
p.stim_cond_duration = zeros(p.stim_cond_count, 1);  % Condition duration in seconds, 0 for unlimited duration.
p.stim_cond_repeat_condition = zeros(p.stim_cond_count, 1); % Go to this condition when duration expires.
p.stim_cond_repetitions = zeros(p.stim_cond_count, 1);      % Maximum number of repetions before continuing to following condition
p.stim_cond_stim_mode = zeros(p.stim_cond_count, 1);        % 0=No stimulation, Otherwise this is the number of active stim channels.
p.stim_isi = 2 * ones(p.stim_cond_count, p.stim_channel_count);  % Minimum Inter-Stimulus Interval
p.stim_level = zeros(p.stim_cond_count, p.stim_channel_count);   % Stimulation level in uA.
p.stim_width = 0.2 * ones(p.stim_cond_count, p.stim_channel_count);  % Phase width in ms.
p.stim_delay = zeros(p.stim_cond_count, p.stim_channel_count);   % Millisecond Stimulus delay from trigger.
p.stim_pulses = zeros(p.stim_cond_count, p.stim_channel_count);  % Pulses in stimulus pulse train. 0 = Don't deliver or mark the stimulus event.
p.stim_level2 = p.stim_level;  % Optional level and width for the second phase (of the bi-phasic pulse)
p.stim_width2 = p.stim_width;
p.stim_cathodic_first = repmat([1 3 5 0 0 0], p.stim_cond_count, 1); % Pin# for Positive voltage on first phase.  Positive current flows out of this pin on the first phase
p.stim_anodic_first = repmat([2 4 6 0 0 0], p.stim_cond_count, 1);   % Pin# for Negative voltage on first phase. 0 
p.stim_limit_count = 200 * ones(p.stim_cond_count, p.stim_channel_count); % Safe stimulation limit. Exceeding this causes a time-out.
p.stim_limit_window = ones(p.stim_cond_count, p.stim_channel_count);      % Safe as long as stimulations <= count inside each window (given in seconds)
p.stim_limit_timeout = 5 * ones(p.stim_cond_count, p.stim_channel_count); % Number of seconds to stop stimulation if safety limit exceeded.

p.stim_req_count = 16; % Maximum allowed stim requirements.
p.stim_req_active_cond = zeros(p.stim_req_count, 1);%0=free, bits for up to 8 conditions (if sent as a byte)
p.stim_req_active_stim = ones(p.stim_req_count, 1); % 0=continuation. Bits for up to 8 stimulus channels (if sent as a byte)
p.stim_req_source = zeros(p.stim_req_count, 1);     % Source variable
p.stim_req_type = zeros(p.stim_req_count, 1);       % Type variable
p.stim_req_n1 = ones(p.stim_req_count, 1);          % Count variable
p.stim_req_t1 = zeros(p.stim_req_count, 1);         % Time variable

% FPGA testing.  This is here so that the last set of FPGA testing
% values will be saved to the settings file.  Normal users should not
% have to see this information in the GUI. Set ui.show_test_tab = 0 at
% the beginning of this function.

p.test_count = 16;
for item = 1:p.test_count
    p.test_name{item} = ['test' num2str(item)];
end
p.test_address = zeros(p.test_count, 1);
p.test_value = zeros(p.test_count, 1);
ui.last_test_index = 1;

% Command flags

ui.command = 0;  % 'g' to get parameters from neurochip, 's' to set paramaters to neurochip, 0 if neither is active.
ui.command_set = 0;     % For pending set parameter commands
ui.command_get = 0;     % For pending get parameter commands
ui.command_record = 0;  % For pending record commands.

% packet header marker byte.

ui.packet_marker = uint8(hex2dec('A5'));
ui.packet_check = uint32(hex2dec('A5A5A5A5'));

% Initialize communication ports.

ui.error = [];  % Normally [] for no error condition.  Set to a string for error.
ui.irport = []; % Serial port object
ui.connection_status = 0; % 0=disconnected, 1=connecton running, 2=disconnect.
ui.irport_status = 1; % 1 = green (last comm succeded, 2 = warning, 3 = error.
ui.status_color = {[0 1 0]; [1 1 0]; [1 0 0]}; % {green; yellow; red}
ui.stream_flag = 0;
ui.stream_time = clock;

% Create figure

ui.figure_position = [2 2 30 20];
ui.units = 'centimeters';
ui.figure = figure('MenuBar','None', 'ToolBar','None', 'Name',mfilename, ...
        'Units',ui.units, 'Position',ui.figure_position, ...
        'NumberTitle','Off', 'Renderer','painters');

% Parameters for layout

ui.dx = 0.25; % Horizontal space between GUI elements.
ui.dy = 0.15; % Vertical space between GUI elements.
ui.button_width = 2;
ui.button_height = 0.75;
ui.text_width = 2.25;
ui.edit_width = 2;
ui.popup_width = 2;
ui.line_height = 0.5;
ui.ty = 0.05;   % Text y adjustment relative to edit box
ui.tx = 0.1;   % Text x adjustment relative to edit box
ui.fontsize = 12;
ui.panel_height = 3.6;
ui.panel_width = 15;
ui.visible = 'on';
ui.test1 = 0;
ui.test2 = 0;

%%%
% Create main control panel

ui.x = ui.dx;
ui.y = ui.figure_position(4) - ui.dx;
ui.control_panel = new_panel('');

ui.y = ui.panel_height - ui.dx;
ui.connect_button = new_button('Connect', @ConnectButton); % Connect/Disconnect
[ui.x, ui.y] = uicoord(ui.connect_button, 'right', 'top');
ui.comport_popup = new_popup({'Com1'; 'Com2'; 'Com3'; 'Com4'; 'Com5'; 'Com6'; 'Com7'; 'Com8'; 'Com9'; 'Com10'; 'Com11'; 'Com12'}, @DefaultCallback);
[ui.x, ui.y] = uicoord(ui.comport_popup, 'right', 'bottom');
ui.y = ui.y - ui.dy;
ui.text_version = new_text('');
ui.x = uicoord(ui.text_version, 'right');
ui.y = uicoord(ui.connect_button, 'top');
ui.set_button = new_button('Set', @SetButton); % Set Parameter button
[ui.x, ui.y] = uicoord(ui.set_button, 'right', 'top');
ui.record_button = new_button('Record', @RecordButton); % Stop/Record
set(ui.record_button, 'Enable', 'Off');
[ui.x, ui.y] = uicoord(ui.set_button, 'left', 'bottom');
ui.y = ui.y - ui.line_height - ui.dy;
ui.text_session_info = new_text('ID: 0');
[ui.x, ui.y] = uicoord(ui.record_button, 'left', 'bottom');
ui.y = ui.y - ui.line_height - ui.dy;
ui.text_record_info = new_text('SD: 0');
[ui.x, ui.y] = uicoord(ui.record_button, 'right', 'top');
ui.edit_width = 1.5;
ui.text_width = 1;
ui.x = ui.x + .4;
ui.uv_scale = new_editbox('100', @UpdateCallback, 'uV');
[ui.x, ui.y] = uicoord(ui.uv_scale, 'left', 'bottom');
ui.x = ui.x - ui.text_width;
ui.ms_scale = new_editbox('20', @UpdateCallback, 'ms');
ui.text_width = 2.25;
ui.edit_width = 2;
[ui.x, ui.y] = uicoord(ui.connect_button, 'left', 'bottom');
ui.subject = new_editbox('Test', @DefaultCallback, 'Subject:');
%[ui.x, ui.y] = uicoord(ui.subject, 'right', 'bottom');
ui.y = ui.y - ui.line_height - ui.dy;
ui.x = ui.x - ui.tx;
ui.settings_text = new_text('Settings File:');
[ui.x, ui.y] = uicoord(ui.settings_text, 'right', 'bottom');
ui.text_width = 12;
ui.x = ui.x - ui.tx;
ui.settings_file = new_text(p.file_name);
set(ui.settings_text, 'HorizontalAlignment', 'right');
ui.text_width = 2.5;
ui.x = ui.dx;
%ui.y = ui.y - ui.line_height;
ui.save_button = new_button('Save...', @SaveButton);     % Save settings button
[ui.x, ui.y] = uicoord(ui.save_button, 'right', 'top');
ui.open_button = new_button('Restore...', @OpenButton);  % Restore settings button
[ui.x, ui.y] = uicoord(ui.open_button, 'right', 'top');
ui.menu_button = new_button('Menu', @MenuBarButton);     % Show menubar button
[ui.x, ui.y] = uicoord(ui.menu_button, 'right', 'top');
ui.savefolder_button = new_button('Download...', @SaveFolderButton);  % Download SDCard file button
[ui.x, ui.y] = uicoord(ui.savefolder_button, 'right', 'top');
ui.stop_button = new_button('Condition Stop', @StopButton);      % Test button for test code.
item = get(ui.stop_button, 'Position');
item(3) = item(3) * 1.5;
item(4) = item(4) * 1.25;
set(ui.stop_button, 'Position', item, 'Enable', 'off');
[ui.x, ui.y] = uicoord(ui.menu_button, 'right', 'top');
ui.x = 12.6;
ui.clear_button = new_button('Clear', @ClearButton);     % Clear sweeps button

% Test button for running TestButton code.
[ui.x, ui.y] = uicoord(ui.clear_button, 'left', 'top');
ui.y = ui.y + ui.button_height + ui.dy;
ui.test_button = new_button('Test', @TestButton);
if (ui.show_test_button == 0)
    set(ui.test_button, 'Visible', 'off'); % Hide test button.
end

%%%
% Create buttons for tabbed panels.

ui.current_panel = ui.figure;
ui.x = 2 * ui.dx;
ui.tab_names = {'Record', 'Event', 'Stim', 'Z Check', 'Test'};
for item = 1:length(ui.tab_names)
    ui.y = uicoord(ui.control_panel, 'bottom') - .1;    
    ui.tab_buttons(item) = new_button(ui.tab_names{item}, @TabSelect);
    if item > 1
        set(ui.tab_buttons(item), 'BackgroundColor', get(ui.figure, 'Color'));
    end
    ui.x = uicoord(ui.tab_buttons(item), 'right');
    if strcmp(ui.tab_names{item}, 'Test') && (ui.show_test_tab == 0)
        set(ui.tab_buttons(item), 'Visible', 'off');
    end
end

%%%
% Create tabbed panels for parameters

ui.panel_height = uicoord(ui.tab_buttons(1), 'bottom') - ui.dx;
ui.x = ui.dx;
for item = 1:length(ui.tab_names)
    ui.y = uicoord(ui.tab_buttons(1), 'bottom');
    ui.tab_panels(item) = new_panel('');
    set(ui.tab_buttons(item), 'Userdata', ui.tab_panels(item));
    ui.visible = 'Off';
end
ui.visible = 'on';
ui.tab_cover = uicontrol('Style','frame', ...
    'Parent', ui.figure, 'Units',ui.units, ...
    'Position',[ui.x, ui.y, ui.button_width, 0.02]);

%%%
% Channel Recording Parameters

ui.x = ui.dx;
ui.y = ui.panel_height - ui.dx;
ui.current_panel = ui.tab_panels(1);
ui.lower_bands = {'500 Hz'; '300 Hz'; '250 Hz'; ...
    '200 Hz'; '150 Hz'; '100 Hz'; '75 Hz'; '50 Hz'; '30 Hz'; '25 Hz'; ...
    '20 Hz'; '15 Hz'; '10 Hz'; '7.5 Hz'; '5.0 Hz'; '3.0 Hz'; '2.5 Hz'; ...
    '2.0 Hz'; '1.5 Hz'; '1.0 Hz'; '0.75 Hz'; '0.50 Hz'; '0.30 Hz';...
    '0.25 Hz'; '0.10 Hz'};
ui.popup_width = 2;
ui.lower_bandwidth_popup = new_popup(ui.lower_bands, @RecordParam, 'Bandwidth:', 'to');
[ui.x, ui.y] = uicoord(ui.lower_bandwidth_popup, 'right', 'top');
ui.x = ui.x + 0.3;
ui.upper_bands = {'20 kHz'; '15 kHz'; '10 kHz'; ...
    '7.5 kHz'; '5.0 kHz'; '3.0 kHz'; '2.5 kHz'; '2.0 kHz'; '1.5 kHz'; ...
    '1.0 kHz'; '750 Hz'; '500 Hz'; '300 Hz'; '250 Hz'; '200 Hz'; ...
    '150 Hz'; '100 Hz'};
ui.upper_bandwidth_popup = new_popup(ui.upper_bands, @RecordParam);

ui.dsp_high_pass_table = [0 -1 0.1103 0.04579 0.02125 0.01027 0.005053 ...
    0.002506 0.001248 0.0006229 0.0003112 0.0001555 0.00007773 0.00003886 ...
    0.00001943 0.000009714 0.000004857];
ui.dsp_high_pass_info{1} = 'Off';
ui.dsp_high_pass_info{2} = 'Differentiator';
for item = 3:8
    ui.dsp_high_pass_info{item} = [num2str(ceil(ui.dsp_high_pass_table(item) * 5000)) ' Hz at 5k, ' num2str(ceil(ui.dsp_high_pass_table(item) * 20000)) ' Hz at 20k'];
end
for item = 8:length(ui.dsp_high_pass_table)
    ui.dsp_high_pass_info{item} = [num2str(ui.dsp_high_pass_table(item) * 5000,'%7.2g') ' Hz at 5k, ' num2str(ui.dsp_high_pass_table(item) * 20000,'%7.2g') ' Hz at 20k'];
end
ui.x = ui.dx;
ui.popup_width = 4.5;
ui.dsp_high_pass_popup = new_popup(ui.dsp_high_pass_info, @RecordParam, 'Offset Removal:', 'DSP high pass');

ui.popup_width = 2;
ui.intan_channels_popup = new_popup({'16'; '32'}, @RecordParam, 'Intan Channels:');
ui.remote_channels_popup = new_popup({'None'; '1'; '2'; '3'; '4'; '5'; '6'; '7'; '8'}, @RecordParam, 'Remote Channels:');
if ui.show_remote_channels == 0
    set(ui.remote_channels_popup, 'visible', 'off');
    set(get(ui.remote_channels_popup, 'UserData'), 'visible', 'off');
    ui.y = ui.y + ui.line_height;
end

ui.x = 1;
ui.y = ui.y - ui.line_height;
uipanel(ui.figure, 'Title', 'Channels', 'Visible', ui.visible, ...
    'FontSize', 9, 'Units', ui.units, 'Parent', ui.current_panel, ...
    'Position', [ui.x - ui.dx, ui.y - 11.5, 13, 11.5]);
ui.y = ui.y - ui.line_height;
ui.chan_list = uicontrol('Style','listbox', 'Tag', 'Channel List',...
    'Max', 3, 'Min', 1', 'Parent', ui.current_panel, 'Units', ui.units, ...
    'Position',[ui.x,ui.y-10.5-ui.dy,4,10.5],...
    'HandleVisibility','callback', 'Callback',@ChanListChanged, ...
    'Background', 'white', 'String',p.channel_names);
ui.edit_width = 2.5;
ui.x = 5;
ui.chan_name = new_editbox('Ch1', @RecordParam, 'Channel Name:');
set(ui.chan_name, 'HorizontalAlignment', 'left');
ui.chan_rate_conversion = [0 5000 10000 20000];
ui.popup_width = 2.5;
ui.chan_rate = new_popup({'Off'; '5 kHz'; '10 kHz'; '20 kHz'}, @RecordParam, 'Sample Rate:');
ui.popup_width = 3.5;
ui.stim_aux3_input = new_popup({'Tissue Ground', 'A1 Voltage Check', 'A0 Voltage Check', ...
    'B1 Voltage Check', 'B0 Voltage Check', 'C1 Voltage Check', 'C0 Voltage Check', ...
    'Tissue Gound (NCal)', 'None'}, @ConditionParam, 'Aux3 Input:');

%%%
% Event Parameters

ui.x = ui.dx;
ui.y = ui.panel_height - ui.dx;
ui.edit_width = 1.5;
ui.current_panel = ui.tab_panels(2);
ui.popup_width = 3;
ui.event_method = new_popup({'Off'; 'Discriminator'; 'Periodic'; 'Discrim + Periodic (on next Event channel)'}, @EventParam, ' ');
ui.popup_width = 2.3;
ui.y = ui.panel_height - ui.dx;
ui.event_list = new_popup({'Event 0'; 'Event 1'; ...
    'Event 2'; 'Event 3'; 'Event 4'; ...
    'Event 5'; 'Event 6'; 'Event 7'}, @EventSelect);
ui.popup_width = 3;
ui.event_source1 = new_popup([{'0'}; p.channel_names], @EventParam, 'Source:', ' -');
[ui.x, ui.y] = uicoord(ui.event_source1, 'right', 'top');
ui.x = ui.x + 0.3;
ui.event_source2 = new_popup([{'0'}; p.channel_names], @EventParam);
ui.x = ui.dx;
ui.event_transform = new_popup({'None'; 'Rectify'}, @EventParam, 'Transform:');
ui.event_lower_bandwidth = new_editbox('10', @EventParam, 'Bandwidth:', 'to');
[ui.x, ui.y] = uicoord(ui.event_lower_bandwidth, 'right', 'top');
ui.x = ui.x + 0.3;
ui.event_upper_bandwidth = new_editbox('200', @EventParam, [], 'Hz (0,0 = off)');

ui.x = ui.dx;
ui.y = ui.y - ui.dx;
ui.event_threshold = new_editbox('0.5', @EventParam, 'Threshold:', 'uV');
ui.event_window1_delay = new_editbox('5', @EventParam, 'Window1 Delay:', 'ms');
ui.event_window1_max = new_editbox('2', @EventParam, 'Max1:', 'uV');
ui.event_window1_min = new_editbox('1', @EventParam, 'Min1:', 'uV');
[ui.x, ui.y] = uicoord(ui.event_window1_delay, 'right', 'top');
ui.x = ui.x + 1;
ui.event_window2_delay = new_editbox('10', @EventParam, 'Window2 Delay:', 'ms');
ui.event_window2_max = new_editbox('-1', @EventParam, 'Max2:', 'uV');
ui.event_window2_min = new_editbox('-2', @EventParam, 'Min2:', 'uV');

ui.x = ui.dx;
ui.y = ui.y - ui.dx;
ui.event_refractory_period = new_editbox('5', @EventParam, 'Refractory Period:', 'ms');
[ui.x, ui.y] = uicoord(ui.event_refractory_period, 'right', 'top');
ui.x = ui.x + 0.25;
item = ui.text_width;
ui.text_width = 1;
ui.event_refractory_count = new_editbox('1', @EventParam, 'after', 'events');
[ui.x, ui.y] = uicoord(ui.event_refractory_count, 'right', 'top');
ui.x = ui.x + 0.75;
ui.text_width = 1.25;
ui.event_refractory_window = new_editbox('5', @EventParam, 'inside of', 'ms');
ui.text_width = item;

ui.x = ui.dx;
ui.y = ui.y - ui.dx - ui.line_height;
new_text('_____________');
ui.x = ui.dx;
ui.y = ui.y - ui.line_height;
new_text('Periodic Events');
ui.x = ui.dx;
ui.event_interval_min = new_editbox('5', @EventParam, 'Minimum Interval:', 'ms');
ui.event_ratio = new_editbox('1', @EventParam, 'Event Ratio 1:');
[ui.x, ui.y] = uicoord(ui.event_interval_min, 'right', 'top');
ui.x = ui.x + 1;
ui.event_interval_max = new_editbox('Off', @EventParam, 'Maximum Interval:', 'ms (0 = off)');
ui.event_interval_exp = new_editbox('Off', @EventParam, 'Poisson Rate:', 'Hz (0 for Uniform)');

ui.event_unit_text = zeros(1,6); % List of text objects denoting scales for Event Sources.
item = get(ui.event_threshold, 'UserData'); % Access to 'uV' units text object.
ui.event_unit_text(1) = item(2);
item = get(ui.event_window1_max, 'UserData'); % Access to 'uV' units text object.
ui.event_unit_text(2) = item(2);
item = get(ui.event_window1_min, 'UserData'); % Access to 'uV' units text object.
ui.event_unit_text(3) = item(2);
item = get(ui.event_window2_max, 'UserData'); % Access to 'uV' units text object.
ui.event_unit_text(4) = item(2);
item = get(ui.event_window2_min, 'UserData'); % Access to 'uV' units text object.
ui.event_unit_text(5) = item(2);
item = get(ui.uv_scale, 'UserData'); % Access to 'uV' units text object.
ui.event_unit_text(6) = item(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stimulation Parameters

% Starting condition and artifact suppression
ui.current_panel = ui.tab_panels(3);
ui.x = ui.dx;
ui.y = ui.panel_height - ui.dx;
ui.popup_width = 3.3;
ui.text_width = 2.5;
ui.stim_power = new_popup({'Off'; '10 Volt Compliance'; '60 Volt Compliance'}, @ConditionParam, 'Stimulator Power:');
ui.start_cond = new_popup({'Condition 1'; 'Condition 2'; 'Condition 3'; ...
    'Condition 4'; 'Condition 5'; 'Condition 6'; 'Condition 7'; ...
    'Condition 8'; 'Condition Stop'}, @ConditionParam, 'Starting Condition:');
[ui.x, ui.y] = uicoord(ui.stim_power, 'right', 'top');
ui.x = uicoord(ui.start_cond, 'right');
ui.text_width = 3.5;
ui.popup_width = 4;
ui.stim_fast_settle = new_popup(p.stim_fast_settle_menu, @ConditionParam, 'Intan Fast Settle:');
ui.x = uicoord(ui.start_cond, 'right');
ui.popup_width = 3;
ui.stim_artifact_time = new_editbox('0', @ConditionParam, 'Event Filter Suppression:', 'ms on each pulse');

% Parameters for each condition
ui.text_width = 2;
[ui.x, ui.y] = uicoord(ui.start_cond, 'left', 'bottom');
ui.y = ui.y - ui.line_height;
ui.x = ui.dx;
new_text('_____________');
ui.x = ui.dx;
ui.cond_list = new_popup({'Condition 1'; 'Condition 2'; 'Condition 3'; 'Condition 4'; ...
    'Condition 5'; 'Condition 6'; 'Condition 7'; 'Condition 8'}, @ConditionSelect, [], 'Settings');
ui.y = ui.y - ui.dy;
ui.cond_duration = new_editbox('unlimited', @ConditionParam, 'Duration:', 'sec');
[ui.x, ui.y] = uicoord(ui.cond_duration, 'right', 'top');
ui.x = ui.x + 1;
ui.cond_repeat_cond = new_popup({'Condition 1'; 'Condition 2'; ...
    'Condition 3'; 'Condition 4'; 'Condition 5'; 'Condition 6'; 'Condition 7'; ...
    'Condition 8'; 'Condition Stop'}, @ConditionParam, 'Proceed to:');
[ui.x, ui.y] = uicoord(ui.cond_repeat_cond, 'right', 'top');

ui.cond_repeat_count = new_editbox('0', @ConditionParam, [], 'repetitions');
ui.x = ui.dx;
ui.popup_width = 2;
ui.cond_stim_mode = new_popup({'None'; '1'; '2'; '3'; '4'; '5'; '6'}, @ConditionParam, 'Channels:');

ui.text_width = 1;
ui.edit_width = 1;
ui.popup_width = 1.1;
[ui.x, ui.y] = uicoord(ui.cond_stim_mode, 'left', 'bottom');
ui.y = ui.y - 2* ui.line_height + ui.dy;
ui.x = ui.dx + ui.text_width;
set(new_text('Delay'), 'HorizontalAlignment', 'center'); ui.x = ui.x - ui.dx/2;
set(new_text('Level'), 'HorizontalAlignment', 'center'); ui.x = ui.x - ui.dx/2;
set(new_text('Width'), 'HorizontalAlignment', 'center'); ui.x = ui.x - ui.dx/2;
set(new_text('Pulse'), 'HorizontalAlignment', 'center'); ui.x = ui.x - ui.dx/2;
set(new_text('ISI'), 'HorizontalAlignment', 'center'); ui.x = ui.x - ui.dx/2;
set(new_text('Limit'), 'HorizontalAlignment', 'center'); ui.x = ui.x - ui.dx/2;
set(new_text('Limit'), 'HorizontalAlignment', 'center'); ui.x = ui.x - ui.dx/2;
set(new_text('Time-'), 'HorizontalAlignment', 'center'); ui.x = ui.x - ui.dx/2;
set(new_text('Stim'), 'HorizontalAlignment', 'center'); ui.x = ui.x - ui.dx/2;
set(new_text('Diff'), 'HorizontalAlignment', 'center'); ui.x = ui.x - ui.dx/2;
ui.visible = 'off';
ui.phase2_uA_text = new_text('2nd uA');
set(ui.phase2_uA_text, 'HorizontalAlignment', 'center'); ui.x = ui.x - ui.dx/2;
ui.phase2_width_text = new_text('Width');
set(ui.phase2_width_text, 'HorizontalAlignment', 'center'); ui.x = ui.x - ui.dx/2;
ui.visible = 'on';
ui.phase2_checkbox = uicontrol('Style','checkbox',...
    'Parent', ui.current_panel, 'Units', ui.units, ...
    'Position',[ui.x - 2 * ui.text_width, ui.y + ui.line_height, ui.button_width, ui.button_height],...
    'HandleVisibility','callback', 'Callback',  @ConditionParam, ...
    'String', '2nd Phase', 'Value', 0);
ui.y = ui.y - ui.line_height + ui.dy;
ui.x = ui.dx + ui.text_width;
set(new_text('ms'), 'HorizontalAlignment', 'center'); ui.x = ui.x - ui.dx/2;
set(new_text('uA'), 'HorizontalAlignment', 'center'); ui.x = ui.x - ui.dx/2;
set(new_text('ms'), 'HorizontalAlignment', 'center'); ui.x = ui.x - ui.dx/2;
set(new_text('count'), 'HorizontalAlignment', 'center'); ui.x = ui.x - ui.dx/2;
set(new_text('ms'), 'HorizontalAlignment', 'center'); ui.x = ui.x - ui.dx/2;
set(new_text('count'), 'HorizontalAlignment', 'center'); ui.x = ui.x - ui.dx/2;
set(new_text('sec'), 'HorizontalAlignment', 'center'); ui.x = ui.x - ui.dx/2;
set(new_text('out sec'), 'HorizontalAlignment', 'center'); ui.x = ui.x - ui.dx/2;
set(new_text('Pin'), 'HorizontalAlignment', 'center'); ui.x = ui.x - ui.dx/2;
set(new_text('Pin'), 'HorizontalAlignment', 'center'); ui.x = ui.x - ui.dx/2;

ui.x = ui.dx;
ui.y = ui.y + ui.dy;
for item = 1:p.stim_channel_count
  ui.stim_delay(item) = new_editbox('0', @ConditionParam, ['Stim ' char(48+item) ':']);
end
[ui.x, ui.y] = uicoord(ui.stim_delay(1), 'right', 'top');
ui.x = ui.x - ui.dx/2;
for item = 1:p.stim_channel_count
    ui.stim_level(item) = new_editbox('0', @ConditionParam);
end
[ui.x, ui.y] = uicoord(ui.stim_level(1), 'right', 'top');
ui.x = ui.x - ui.dx/2;
for item = 1:p.stim_channel_count
    ui.stim_width(item) = new_editbox('0.2', @ConditionParam);
end
[ui.x, ui.y] = uicoord(ui.stim_width(1), 'right', 'top');
ui.x = ui.x - ui.dx/2;
for item = 1:p.stim_channel_count
    ui.stim_pulses(item) = new_editbox('0', @ConditionParam);
end
[ui.x, ui.y] = uicoord(ui.stim_pulses(1), 'right', 'top');
ui.x = ui.x - ui.dx/2;
for item = 1:p.stim_channel_count
    ui.stim_isi(item) = new_editbox('1', @ConditionParam);
end
[ui.x, ui.y] = uicoord(ui.stim_isi(1), 'right', 'top');
ui.x = ui.x - ui.dx/2;
for item = 1:p.stim_channel_count
    ui.stim_limit_count(item) = new_editbox('200', @ConditionParam);
end
[ui.x, ui.y] = uicoord(ui.stim_limit_count(1), 'right', 'top');
ui.x = ui.x - ui.dx/2;
for item = 1:p.stim_channel_count
    ui.stim_limit_window(item) = new_editbox('1', @ConditionParam);
end
[ui.x, ui.y] = uicoord(ui.stim_limit_window(1), 'right', 'top');
ui.x = ui.x - ui.dx/2;
for item = 1:p.stim_channel_count
    ui.stim_limit_timeout(item) = new_editbox('5', @ConditionParam);
end

% Cathode/Anode selection
[ui.x, ui.y] = uicoord(ui.stim_limit_timeout(1), 'right', 'top');
ui.x = ui.x - ui.dx/2;
for item = 1:p.stim_channel_count
    ui.stim_cathode(item) =  new_popup({'Off'; 'A1'; 'A0'; 'B1'; 'B0'; 'C1'; 'C0'}, @ConditionParam);
    if item <= 3
        set(ui.stim_cathode(item), 'Value', item * 2);
    end
end
[ui.x, ui.y] = uicoord(ui.stim_cathode(1), 'right', 'top');
ui.x = ui.x - ui.dx/2;
for item = 1:p.stim_channel_count
    ui.stim_anode(item) =  new_popup({'Off'; 'A1'; 'A0'; 'B1'; 'B0'; 'C1'; 'C0'}, @ConditionParam);
    if item <= 3
        set(ui.stim_anode(item), 'Value', item * 2 + 1);
    end
end

% Optional 2nd Phase uA Level and ms width.
ui.visible = 'off';
[ui.x, ui.y] = uicoord(ui.stim_anode(1), 'right', 'top');
ui.x = ui.x - ui.dx/2;
for item = 1:p.stim_channel_count
    ui.stim_level2(item) = new_editbox('0', @ConditionParam);
end
[ui.x, ui.y] = uicoord(ui.stim_level2(1), 'right', 'top');
ui.x = ui.x - ui.dx/2;
for item = 1:p.stim_channel_count
    ui.stim_width2(item) = new_editbox('0.2', @ConditionParam);
end
ui.visible = 'on';

% Stimulation requirements
ui.x = ui.dx;
ui.y = ui.y - ui.line_height;
uipanel(ui.figure, 'Title', 'Stimulation Requirements', 'Visible', ui.visible, ...
    'FontSize', 9, 'Units', ui.units, 'Parent', ui.current_panel, ...
    'Position', [ui.x, ui.y - 5.25, 14.25, 5.25]);

ui.y = ui.y - ui.line_height;
ui.x = ui.x + ui.dx;
ui.req_list = uicontrol('Style','listbox', 'Tag', 'Requirement List',...
    'Max', 1, 'Min', 1', 'Parent', ui.current_panel, 'Units', ui.units, ...
    'Position',[ui.x,ui.y-3-ui.dy,13.5,3],...
    'HandleVisibility','callback', 'Callback',@ReqSelect, ...
    'Background', 'white', 'String', repmat({'1) Off'}, p.stim_req_count, 1));

ui.y = ui.y - 3 - ui.dy;
ui.popup_width = 3;
ui.text_width = 1.1;
ui.req_cond = new_popup({'None'; '  Condition 1'; '  Condition 2'; '  Condition 3'; ...
    '  Condition 4'; '  Condition 5'; '  Condition 6'; '  Condition 7'; '  Condition 8'}, @ReqParam, 'In' );

[ui.x, ui.y] = uicoord(ui.req_cond, 'right', 'top');
ui.text_width = 1;
[ui.x, ui.y] = uicoord(ui.req_cond, 'right', 'top');
ui.popup_width = 2;
ui.req_source = new_popup({'Event 1'; 'Event 2'; 'Event 3'; 'Event 4'; ...
    'Event 5'; 'Event 6'; 'Event 7'; 'Event 8'}, @ReqParam, 'When');

[ui.x, ui.y] = uicoord(ui.req_source, 'right', 'top');
ui.popup_width = 3.75;
ui.req_type = new_popup({'occurs'; 'count > N'; 'count <= N'; 'Source inside Windows'; 'Source outside Windows'}, @ReqParam);

ui.text_width = 1.1;
ui.x = 2 * ui.dx;
ui.popup_width = 3;
ui.req_output_chans = new_popup({'None'; '  Stim 1'; '  Stim 2'; ...
    '  Stim 3'; '  Stim 4'; '  Stim 5'; '  Stim 6'}, @ReqParam, 'Trigger');

ui.text_width = 1.1;
[ui.x, ui.y] = uicoord(ui.req_type, 'left', 'bottom');
ui.req_t1 = new_editbox('0', @ReqParam, 'for last', 'ms');
ui.text_width = 1;
[ui.x, ui.y] = uicoord(ui.req_type, 'right', 'top');
ui.req_n1 = new_editbox('0', @ReqParam, 'N =');

%%%
% Impedance (Z Check) panel

ui.x = ui.dx;
ui.y = ui.panel_height - ui.dx;
ui.current_panel = ui.tab_panels(4);

ui.text_width = 3;
ui.popup_width = 3;
if (ui.show_10MOhm_scale == 0)
    % Normal impedance test popup menu.
    ui.z_check_menu = new_popup({'Off'; '100 kOhm scale'; '1 MOhm scale'}, @ZCheckParam, 'Impedance Check:');
else
    % 10 MOhm scale seems off by a factor of 2x to 3x, but it can be
    % enabled by setting ui.show_10MOhm_scale to 1 at the top of this file.
    ui.z_check_menu = new_popup({'Off'; '100 kOhm scale'; '1 MOhm scale'; '10 MOhm scale'}, @ZCheckParam, 'Impedance Check:');
end
ui.z_check_output = uicontrol('Style','edit', 'Tag', 'ZCheck Output', 'HorizontalAlignment', 'left', ...
    'Max', 3, 'Min', 1', 'Parent', ui.current_panel, 'Units', ui.units, ...
    'Position',[ui.x,ui.y-13.25-ui.line_height,10,13.5], 'Background', 'white', 'String', '');

%%%
% Test panel

ui.x = ui.dx;
ui.y = ui.panel_height - ui.dx - ui.line_height;
ui.current_panel = ui.tab_panels(5);
ui.edit_width = 3;
for item = 1:p.test_count
    ui.test_name(item) = new_editbox(['test' num2str(item)], @TestParam);
end
[ui.x, ui.y] = uicoord(ui.test_name(1), 'right', 'top');
ui.edit_width = 1.5;
for item = 1:p.test_count
    ui.test_address(item) = new_editbox('0x0000', @TestParam);
end
ui.edit_width = 2.25;
[ui.x, ui.y] = uicoord(ui.test_address(1), 'right', 'top');
for item = 1:p.test_count
    ui.test_value(item) = new_editbox('0x00000000', @TestParam);
end
[ui.x, ui.y] = uicoord(ui.test_value(1), 'right', 'top');
ui.button_height = ui.line_height;
ui.button_width = 1.5;
for item = 1:p.test_count
    ui.test_fpga_button(item) = new_button('Write', @TestFPGAButton);
end
[ui.x, ui.y] = uicoord(ui.test_fpga_button(1), 'right', 'top');
ui.button_height = ui.line_height;
ui.button_width = 1.5;
for item = 1:p.test_count
    ui.test_fpga_button(item+p.test_count) = new_button('Read', @TestFPGAButton);
end
[ui.x, ui.y] = uicoord(ui.test_fpga_button(p.test_count+1), 'right', 'top');
for item = 1:p.test_count
    ui.test_response(item) = new_editbox('0x00000000', @DefaultCallback);
end
[ui.x, ui.y] = uicoord(ui.test_name(1), 'left', 'top');
ui.y = ui.y - ui.tx;
new_text('Name');
[ui.x, ui.y] = uicoord(ui.test_address(1), 'left', 'top');
ui.y = ui.y - ui.tx;
new_text('Address');
[ui.x, ui.y] = uicoord(ui.test_value(1), 'left', 'top');
ui.y = ui.y - ui.tx;
new_text('Value');
[ui.x, ui.y] = uicoord(ui.test_response(1), 'left', 'top');
ui.y = ui.y - ui.tx;
new_text('Response');

%%%
% Create Plot Panel

ui.visible = 'on';
ui.plot_panel = uipanel(ui.figure, 'Units', ui.units, ...
	'Position',[ui.panel_width + 2*ui.dx, ui.dx, ui.figure_position(3) - ui.panel_width - 3*ui.dx, ui.figure_position(4) - 2*ui.dx]);
ui.plot = subplot('Position',[0.1 0.71 .85 .28]);
if isprop(ui.plot,'SortMethod')
    set(ui.plot, 'Parent', ui.plot_panel, 'SortMethod','childorder');
else
    set(ui.plot, 'Parent', ui.plot_panel, 'DrawMode', 'fast');
end
for item = 1:p.event_count
    ui.chplot(item) = subplot('Position', [0.1, .7 - item/12, .85, .06]);
    if isprop(ui.plot,'SortMethod')
        set(ui.chplot(item), 'Parent', ui.plot_panel, 'XTickMode', 'manual', 'SortMethod','childorder');
    else
        set(ui.chplot(item), 'Parent', ui.plot_panel, 'XTickMode', 'manual', 'DrawMode', 'fast');
    end
    ylabel(['E' num2str(item-1)], 'Rotation', 0);
    if item == p.event_count
        set(ui.chplot(item), 'XTickMode', 'auto');
    end
end

% Storage from sweeps shown in the first plot.
ui.parameter_set_status = 0;  % Set to 1 if neurochip status block indicates a successfull parameter_set command.
ui.status_version = 0;
ui.status_msecs = 0;
ui.status_summaries = 0;
ui.status_str = '';
ui.last_sweep_index = 0;
ui.current_sweep = zeros(1, 500);
ui.sweep_info = 256 + 1;  % Last argument sent with 'i' command.
ui.max_sweeps = 20;
ui.max_events = 8;  
ui.sweeps = zeros(ui.max_sweeps, 500);
ui.sweep_ok = zeros(ui.max_sweeps, 1);   % True if sweep(i) data is valid

% Circluar buffer for the Event summaries sent by the neurochip.
ui.event_summary_index = 1;  % last summary index (the one that may still be in progress)
ui.event_summary_count = 100;
ui.event_summary_min = zeros(ui.event_summary_count,p.event_count);
ui.event_summary_max = zeros(ui.event_summary_count,p.event_count);
ui.event_summary_event = zeros(ui.event_summary_count,p.event_count);
ui.event_summary_stim = zeros(ui.event_summary_count,p.stim_channel_count);

% Timer for sending information to the neurochip

ui.command_timer = timerfind('Tag', 'NCUICommandTimer');
if isempty(ui.command_timer)
    ui.command_timer = timer('TimerFcn', @CommandTimerFcn, 'ExecutionMode', 'fixedDelay', 'Period' ,1, 'Tag', 'NCUICommandTimer');
else
    set(ui.command_timer, 'TimerFcn', @CommandTimerFcn, 'ExecutionMode', 'fixedDelay', 'Period' ,1, 'Tag', 'NCUICommandTimer');
end

%%%
% Initialize

tic;
load_last_settings();
update_all_controls();
TabSelect(ui.tab_buttons(1), 0);
irport_init();
update_plots();

% At this point control returns to Matlab and all further actions
% are handled through callbacks from the GUI
%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nested functions follow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function varargout = uicoord(handle, varargin)
        % Returns coordinates for the edges of the space allocated
        % for a GUI element.
        position = get(handle, 'Position');
        varargout = cell(nargin-1,1);
        for i = 1:nargin-1
            switch varargin{i}
                case 'top'
                    varargout{i} = position(2) + position(4) + ui.dy;
                case 'bottom'
                    varargout{i} = position(2);
                case 'left'
                    varargout{i} = position(1);
                case 'right'
                    varargout{i} = position(1) + position(3) + ui.dx;
                otherwise
                    varargout{i} = 0;
            end
        end
    end

    function panel = new_panel(panel_name)
        % Add a new panel to the GUI. Top = ui.y, Left = ui.x
        ui.current_panel = uipanel(ui.figure, 'Title', panel_name, ...
            'FontSize', ui.fontsize, 'Units', ui.units, ...
            'Visible', ui.visible, ...
            'Position', [ui.x, ui.y - ui.panel_height, ui.panel_width, ui.panel_height]);
        panel = ui.current_panel;
    end

	function button = new_button(button_name, callback)
        % Add button to GUI below last element.
        % Leaves ui.x,ui.y at left,bottom of the created button.
        ui.y = ui.y - ui.dy - ui.button_height;
        button = uicontrol('Style','pushbutton',...
            'Parent', ui.current_panel, 'Units', ui.units, ...
            'Position',[ui.x, ui.y, ui.button_width, ui.button_height],...
            'HandleVisibility','callback', 'Callback', callback, ...
            'String', button_name, 'Visible', ui.visible);
    end

    function popup = new_popup(cell_names, callback, varargin)
        % Add pop-up menu below last element
        % Leaves ui.x,ui.y at left,bottom of the created pop-up.
        ui.y = ui.y - ui.dy - ui.line_height;
        handles = [];
        dx = 0;
        if nargin > 2
            if ~isempty(varargin{1})
                dx = ui.text_width;
                handles = uicontrol('Style','text', ...
                    'Parent', ui.current_panel, 'Units', ui.units, ...
                    'Position', [ui.x, ui.y, dx - ui.tx, ui.line_height - ui.tx], ...
                    'HorizontalAlignment', 'right', 'Max', 1, 'Min', 1', ...
                    'String', varargin{1}, 'Visible', ui.visible);
            end
            if nargin > 3
                handles(end+1) = uicontrol('Style','text', ...
                    'Parent', ui.current_panel, 'Units', ui.units, ...
                    'Position', [ui.x + ui.tx + dx + ui.popup_width, ui.y, ui.text_width - ui.tx, ui.line_height - ui.tx],...
                    'HorizontalAlignment', 'left', 'Max', 1, 'Min', 1', ...
                    'String', varargin{2}, 'Visible', ui.visible);
            end
        end
        popup = uicontrol('Style','popupmenu', 'Visible', ui.visible,...
             'Parent', ui.current_panel, 'Units',ui.units, ...
             'Position',[ui.x + dx, ui.y, ui.popup_width, ui.line_height],...
             'HandleVisibility','callback', 'Callback',callback, ...
             'Background', 'white', 'String',cell_names, ...
             'UserData', handles, 'Visible', ui.visible, 'Value', 1);      
    end

    function text = new_text(text_string)
        % Add an uneditable text string at ui.x,ui.y
        % Leaves ui.x,ui.y at the next tab spot.
        text = uicontrol('Style','text', ...
            'Parent', ui.current_panel, 'Units', ui.units, ...
            'Position', [ui.x, ui.y, ui.text_width, ui.line_height], ...
            'HorizontalAlignment', 'left', 'Max', 1, 'Min', 1', ...
            'String', text_string, 'Visible', ui.visible);
        ui.x = ui.x + ui.dx + ui.text_width;
    end

    function editbox = new_editbox(start_text, callback, varargin)
        % Add and edit box below last element
        % Leaves ui.x,ui.y at left,bottom of the created pop-up.
        ui.y = ui.y - ui.dy - ui.line_height;
        handles = [];
        dx = 0;
        if nargin > 2
            if ~isempty(varargin{1})
                dx = ui.text_width;
                handles = uicontrol('Style','text', ...
                    'Parent', ui.current_panel, 'Units', ui.units, ...
                    'Position', [ui.x, ui.y, dx - ui.tx, ui.line_height - ui.ty], ...
                    'HorizontalAlignment', 'right', 'Max', 1, 'Min', 1', ...
                    'String', varargin{1}, 'Visible', ui.visible);
            end
            if nargin > 3
                handles(end+1) = uicontrol('Style','text', ...
                    'Parent', ui.current_panel, 'Units', ui.units, ...
                    'Position', [ui.x + ui.tx + dx + ui.edit_width, ui.y, ui.text_width - ui.tx, ui.line_height - ui.ty],...
                    'HorizontalAlignment', 'left', 'Max', 1, 'Min', 1', ...
                    'String', varargin{2}, 'Visible', ui.visible);
            end
        end
        editbox = uicontrol('Style','edit', 'Visible', ui.visible,...
             'Parent', ui.current_panel, 'Units',ui.units, ...
             'Position',[ui.x + dx, ui.y, ui.edit_width, ui.line_height],...
             'HandleVisibility',' callback', 'Callback', callback, ...
             'Background', 'white', 'String', start_text, ...
             'UserData', handles, 'Visible', ui.visible);
    end        


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions for Plotting Event windows
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function update_plots()
        subplot(ui.plot);
        cla(ui.plot);
        hold(ui.plot, 'on');
        wcolor = [0 0.45 0.75]; % Color for discrimnator windows.
        
        index = get(ui.event_list, 'Value');
        yThreshold = p.event_threshold(index);
        delay1 = p.event_window1_delay(index);
        max1 = p.event_window1_max(index);
        min1 = p.event_window1_min(index);
        delay2 = p.event_window2_delay(index);
        max2 = p.event_window2_max(index);
        min2 = p.event_window2_min(index);
        conversion =  ui.adu2uv; 
        
        source = p.event_source1(index);    
        if source == 0
            source = p.event_source2(index);
        end
        if ismember(source, ui.channel_Accel_index)
            % Handle accelermeter plot separate
            sample1 = round(delay1 / 10) + 400;
            sample2 = round(delay2 / 10) + 400;
            sample1 = max(301, min(500, sample1));
            sample2 = max(301, min(500, sample2));
            xmin = -1000;
            xmax = 1000;
            yscale = getnum(ui.uv_scale);
            xscale = getnum(ui.ms_scale);
            t = (-100:99) * 10;  % Accelerometer samples are 10 ms appart.
    
            conversion = 1.0 / ui.accel2adu; %  Convert to 1000 units = 1g
            if source == ui.channel_Accel_index(5)
                conversion = 1000 / 256;
                sweeps = conversion * sqrt(abs(ui.sweeps * 8));  % Accelerometer Mag channel is ((x*x)+(y*y)+(z*z))/8. where x,y,z are 10-bit values of Axis-X,Y,Z
            elseif source == ui.channel_Accel_index(4)
                conversion = 1.0 / ui.accelTemp2adu;  % Temperature channel has its own conversion to Celcius
                sweeps = conversion * ui.sweeps;
            else
                sweeps = conversion * ui.sweeps;
            end
            
            v1 = sweeps(:, sample1);
            v2 = sweeps(:, sample2);
            sweep_flags = (v1 >= min1) & (v1 < max1) & (v2 >= min2) & (v2 < max2);
            sweep_index = find((sweep_flags == 0) & ui.sweep_ok);
            if ~isempty(sweep_index)
                plot(ui.plot, t, sweeps(sweep_index, end-199:end), 'Color', [.5 .5 .5]);
            end
            sweep_index = find(sweep_flags & ui.sweep_ok);
            if ~isempty(sweep_index)
                plot(ui.plot, t, sweeps(sweep_index, end-199:end), 'Color', [0 0.85 0]);
            end
            line([0, 0], [-yscale, yscale], 'Color', wcolor, 'Parent', ui.plot);
            if yThreshold < 6000
                line([xmin, xmax], [yThreshold, yThreshold], 'Color', wcolor, 'Parent', ui.plot);
            end
            line([delay1, delay1, nan, delay2, delay2], [max1, min1, nan, max2, min2], 'LineWidth', 3, 'Color', wcolor, 'Parent', ui.plot);
            ylabel(ui.plot, ['F' num2str(index)], 'Rotation', 0);
            if xscale > 0
                xlim(ui.plot, [-xscale, xscale]);
            else
                xlim(ui.plot, [xmin, xmax]);
            end
            if yscale > 0
                ylim(ui.plot, [-yscale yscale]);
            end
        else
            ms2samples = ui.ms2fs;  % Conversion factor from ms to number of samples
            if (source > 0)
                ms2samples = p.channel_rate(source) / 1000;
            end
            decimate = double(max(1, bitshift(ui.sweep_info, -8))); % sweep_decimation();
            sample1 = round(delay1 * ms2samples / decimate) + 100;
            sample2 = round(delay2 * ms2samples / decimate) + 100;
            sample1 = max(1, min(500, sample1));
            sample2 = max(1, min(500, sample2));
            xmin = decimate * -5;
            xscale = getnum(ui.ms_scale);
            if xscale > 0
                xmax = xscale;
            else
                xmax = decimate * 20;
            end
            if -xmin > xmax
                xmin = -xmax;
            end
            yscale = getnum(ui.uv_scale);
            t = (-100:399) * decimate / ms2samples;
            
            v1 = ui.sweeps(:, sample1) * conversion;
            v2 = ui.sweeps(:, sample2) * conversion;
            sweep_flags = (v1 >= min1) & (v1 < max1) & (v2 >= min2) & (v2 < max2);
            sweep_index = find((sweep_flags == 0) & ui.sweep_ok);
            if ~isempty(sweep_index)
                plot(ui.plot, t, conversion * ui.sweeps(sweep_index, :), 'Color', [.5 .5 .5]);
            end
            sweep_index = find(sweep_flags & ui.sweep_ok);
            if ~isempty(sweep_index)
                plot(ui.plot, t, conversion * ui.sweeps(sweep_index, :), 'Color', [0 0.85 0]);
            end
            line([0, 0], [-yscale, yscale], 'Color', wcolor, 'Parent', ui.plot);
            if yThreshold < 6000
                line([xmin, xmax], [yThreshold, yThreshold], 'Color', wcolor, 'Parent', ui.plot);
            end
            line([delay1, delay1, nan, delay2, delay2], [max1, min1, nan, max2, min2], 'LineWidth', 3, 'Color', wcolor, 'Parent', ui.plot);
            ylabel(ui.plot, ['F' num2str(index)], 'Rotation', 0);
            xlim(ui.plot, [xmin, xmax]);
            if yscale > 0
                ylim(ui.plot, [-yscale yscale]);
            else
                ylim(ui.plot, [-100,100]);
            end
        end
        
        % Draw each Event plot.
        
        index1 = ui.event_summary_index+1:ui.event_summary_count;
        index2 = 1:ui.event_summary_index;
        t = (-ui.event_summary_count+1:0)' * 0.1;
        
        for i=1:p.event_count
            ax = ui.chplot(i); % The axis in which to plot
            subplot(ax);
            conversion = ui.adu2uv; % Adjust for filter gain
            source = p.event_source1(i);    
            if source == 0
                source = p.event_source2(i);
            end
            if ismember(source, ui.channel_Accel_index)
                if source == ui.channel_Accel_index(4)
                    conversion = 1.0 / ui.accelTemp2adu;  % Temperature has its own conversion
                else
                    conversion = 1.0/ui.accel2adu; %  Convert to 1000 units = 1g
                end
                if source == ui.channel_Accel_index(5)
                    % Magnitude channel is squared and divided by 4
                    conversion = 1000 / 256;
                    dmin = conversion * [sqrt(abs(ui.event_summary_min(index1, i)) * 8); sqrt(abs(ui.event_summary_min(index2, i)) * 8)];
                    dmax = conversion * [sqrt(abs(ui.event_summary_max(index1, i)) * 8); sqrt(abs(ui.event_summary_max(index2, i)) * 8)];
                else
                    dmin = conversion * [ui.event_summary_min(index1, i); ui.event_summary_min(index2, i)];
                    dmax = conversion * [ui.event_summary_max(index1, i); ui.event_summary_max(index2, i)];
                end
            else
                dmin = conversion * [ui.event_summary_min(index1, i); ui.event_summary_min(index2, i)];
                dmax = conversion * [ui.event_summary_max(index1, i); ui.event_summary_max(index2, i)];
            end
            devt = [ui.event_summary_event(index1, i); ui.event_summary_event(index2, i)];
            cla(ax);
            hold(ax, 'on');
            yThreshold = p.event_threshold(i);
            if yThreshold < 6000
                line([-100 0], [yThreshold yThreshold], 'Color', wcolor, 'Parent', ax);
            end
            if sum(dmax - dmin) > 0
                plot(ax, [t, t], [dmin, dmax], 'k'); 
            end
            if sum(devt) > 0
                yl = ylim(ax);
                plot(ax, t, ((yl(2) - yl(1)) / max(devt)) * devt + yl(1), 'b');  % Max rate scaled to max Y-scale with 0Hz = min Y-scale
                h = text(-9.6, yl(1) + 0.8 * (yl(2)-yl(1)), [num2str(10 * sum(devt) / length(devt), 4) ' Hz']);
                set(h, 'Color', 'b', 'Parent', ax);
            end
            if i <= p.stim_channel_count
                dstim = [ui.event_summary_stim(index1, i); ui.event_summary_stim(index2, i)];
                if sum(dstim) > 0
                    yl = ylim(ax);
                    plot(ax, t, ((yl(2) - yl(1)) / max(dstim)) * dstim + yl(1), 'r'); % Max stim rate scaled max Y-scale with 0Hz = min Y-scale
                    h = text(-9.6, yl(1) + 0.25 * (yl(2)-yl(1)), [num2str(10 * sum(dstim) / length(dstim), 4) ' Hz']);
                    set(h, 'Color', 'r', 'Parent', ax);
                end
            end
            if i < p.event_count
                set(ax, 'YLimMode', 'auto', 'XTickLabel', []);
            else
                set(ax, 'YLimMode', 'auto');
            end
            xlim(ax, [-10 0]);
            ylabel(ax, ['E' num2str(i)], 'Rotation', 0);
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions for dealing with the IRPort
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initialize the COM port.
    % Called when initiating a connection to the neurochip Ir port.
    function irport_init()            
        try
            fclose(instrfind); % Try to close all serial port objects
        catch
        end
        comportstrs = get(ui.comport_popup, 'String');
        comport = comportstrs{get(ui.comport_popup, 'Value')};
        ui.irport = instrfind('Tag', 'NC3IRSettings');
        if isempty(ui.irport)
            ui.irport = serial(comport);
        else
            ui.irport = ui.irport(1); % Since instrfind can return multiple ports
        end
        if (isempty(ui.irport))
            errordlg('Unable to create serial port object');
        else
            irport_close();
            set(ui.irport, 'Port', comport, 'Tag','NC3IRSettings', 'BaudRate',57600, ...
                'InputBufferSize',8192, 'OutputBufferSize',8192, 'TimeOut', 1, 'StopBits',1, ...
                'Parity','none', 'DataBits',8, 'FlowControl', 'none');
            irport_open();
        end
    end

    %%%%%%%%%%
    % Opens serial port for IR connection.
    function irport_open
        try
            fopen(ui.irport); % Open IR Port if it is closed.
        catch
            ui.irport = [];
            errordlg('Could not open IrDA port');
        end
    end

    %%%%%%%%%%
    % Closes the IR Port.
    function irport_close
        if ~isempty(ui.irport)
            try
                fclose(ui.irport);
            catch
            end
        end
    end

    %%%%%%%%%%
    % Returns the current communication port in an open state.
    % Returns empty if irport_init has not been called successfully.
    function port = irport()
        if ~isempty(ui.irport)
            if strcmp(get(ui.irport, 'Status'), 'closed')
                irport_open(); % Open IR Port if it is closed.
            end
        end
        port = ui.irport;
    end

    %%%%%%%%%%
    % Returns the value of a ui control String as a number.
    % Performs range coersion if varargin contains [min max].
    function n = getnum(hHandle, varargin)
        style = get(hHandle, 'Style');
        if strcmp(style, 'popupmenu')
            n = get(hHandle, 'Value') - 1; % Popup menues return 0...n-1
        elseif strcmp(style, 'checkbox')
            n = get(hHandle, 'Value');     % Checkboxes return 0=off, 1=on.
        else
            s = get(hHandle, 'String');
            try
                l = length(s);
                if (l>2) && (s(2) == 'x')
                    n = hex2dec(s(3:end));
                else
                    n = str2double(s);
                end
            catch
                n = 0;
            end
            if nargin > 1
                range = varargin{1};
            else
                range = [];
            end
            if (length(range) >= 2)
                if isnan(n) || isempty(n)
                    if (range(1) < 0) && (range(2) > 0)
                        n = 0;
                    elseif range(1) >= 0
                        n = range(1);
                    else
                        n = range(2);
                    end
                    set(hHandle, 'String', num2str(n));
                elseif n < range(1)
                    n = range(1);
                    set(hHandle, 'String', num2str(n));
                elseif n > range(2)
                    n = range(2);
                    set(hHandle, 'String', num2str(n));
                end
            else
                if isempty(n) || isnan(n)
                    n = 0;
                    set(hHandle, 'String', '0');
                end
            end
        end
    end

	%%%%%%%%%%
    % Sets the value of a ui control String as a number.
    % Does nothing if value is empty or NaN.
    % Does range coersion if varargin = [min max].
    function n = setnum(hHandle, value, varargin)
        n = value;
        style = get(hHandle, 'Style');
        if strcmp(style, 'popupmenu')
            n = max(1, floor(n + 1));
            if (n > length(get(hHandle, 'String')))
                n = 1;
            end
            set(hHandle, 'Value', n); % Popup menues transform 0...n-1 to 1..n
        elseif strcmp(style, 'checkbox')
            set(hHandle, 'Value', floor(n)~=0);     % Checkboxes set 0=off, 1=on.
        elseif ~isempty(value)
            if ~isnan(value)
                if nargin > 2
                    range = varargin{1};
                    if (length(range) >= 2)
                        if n < range(1)
                            n = range(1);
                        elseif n > range(2)
                            n = range(2);
                        end
                    end
                end
                set(hHandle, 'String', num2str(n));
            end
        end
    end

    %%%
    % Converts the strings in the bandwidth menu to values in Hertz
    function hertz = str2Hz(string)
        hertz = str2double(string(1:strfind(string, ' ')-1));
        if strfind(string, 'k')
            hertz = hertz * 1000; % For values given in kHz.
        end
    end

    %%%
    % Writes a block of values to the irport.
    % Values are formatted into a stream of unsigned characters based on
    % their class (type).  Values must be either char or integer arrays.
    % Does not check for proper byte alignment required by the ARM CPU.
    % Using little endian byte order for both Intel and SAM4 ARM Cortex4.
    function ir_write(varargin)
        bytes = uint8(zeros(1024,1));
        nbytes = 0;
        for i=1:nargin
            cl = class(varargin{i});
            switch cl
                case {'char', 'int8', 'uint8'}
                    v = double(varargin{i}) + 256;
                    for j=1:length(v)
                        nbytes = nbytes + 1;
                        bytes(nbytes) = uint8(bitand(v(j), 255));
                    end
                case {'int16', 'uint16'}
                    v = double(varargin{i}) + 65536;
                    for j=1:length(v)
                        bytes(nbytes + 1) = uint8(bitand(bitshift(v(j),-8), 255));
                        bytes(nbytes + 2) = uint8(bitand(v(j), 255));
                        nbytes = nbytes + 2;
                    end
                case {'int32', 'uint32'}
                    v = double(varargin{i}) + 4294967296;
                    for j=1:length(v)
                        bytes(nbytes + 1) = uint8(bitand(bitshift(v(j),-24), 255));
                        bytes(nbytes + 2) = uint8(bitand(bitshift(v(j),-16), 255));
                        bytes(nbytes + 3) = uint8(bitand(bitshift(v(j),-8), 255));
                        bytes(nbytes + 4) = uint8(bitand(v(j), 255));
                        nbytes = nbytes + 4;
                    end
                otherwise
                    disp(['Error: ir_write cannot format class: ' cl]);
            end
        end
        if (nbytes > 0)
            try
                port = irport();
                % Write at most 1024 bytes at a time.
                for startbyte = 1:1024:nbytes
                    endbyte = min(startbyte+1023, nbytes);
                    while strcmp(port.TransferStatus, 'write')
                        pause(0.01);
                    end
                    fwrite(port, bytes(startbyte:endbyte), 'uint8', 'async');
                end              
            catch
            end
        end
    end

    %%%
    % Functions for reading/writing data blocks from ir port 
    
    function GetStatusBlock()
        % Read status block recived in ui.r
        % Reference:     
        % ui.event_summary_index = 0;  % last summary index (the one that may still be in progress)       
        % ui.event_summary_count = 100;
        % ui.event_summary_min = zeros(ui.event_summary_count,8);
        % ui.event_summary_max = zeros(ui.event_summary_count,8);
        % ui.event_summary_event = zeros(ui.event_summary_count,8);
        % ui.event_summary_stim = zeros(ui.event_summary_count,6);
        % 	uint8_t header2[4];    // Header for info block  [0] = 'i', [1..3] = PACKET_MARKER
        % 	uint8_t status[4];     // [0] = parameter received ok flag, [1] = current discriminator index, [2] = current condition index.
        % 	uint32_t version;      // Firmware version info
        % 	uint32_t msecs;        // Number of milliseconds since startup.
        % 	uint32_t summaries;    // Number of summaries (last summary may be incomplete).  Low 4 bits = summary index
        % 	int16_t min[16*8];   // offset 20  // 16 summaries at 10 summaries per second = 100ms summary bins.
        % 	int16_t max[16*8];   // offset 276  // Min and max filtered event sample values for each summary.
        % 	uint8_t events[16*8];// offset 532  // Event counts for each summary.
        % 	uint8_t stims[16*6]; // offset 660  // Stimulus counts for each summary.
        %   uint32_t zcheck[32]; // offset 756  // Impedance check values
        % 	uint32_t pad[31];      // space left over (padding to make a 1024 byte packet)
        %   uint16_t session_id; // copy of the session id
        %   uint16_t sweepinfo;   // includes filter index and sweep decimation.
        % 	uint32_t testval1;    // Checksum for summary block
        % 	uint32_t testval2;    // Checksum for summary block
        % 	uint32_t checksum2;    // Checksum for summary block

        ui.parameter_set_status = ui.r(5);
        ui.status_version = ui.r(9:12);
        ui.status_msecs = GetULongReg(12);
        ui.status_summaries = GetULongReg(16);
        sid = GetUWordReg(1008);  % copy of session id
        if (sid > ui.session_id)
            ui.session_id = sid;
        end
        sweepinfo = GetUWordReg(1010);  % filter index and decimation of current filter.
        if (ui.sweep_info ~= sweepinfo)
            ui.sweep_info = sweepinfo;
            ui.sweep_ok(:) = 0;  % Clear sweeps if scale box object was changed.
        end
        icond = ui.r(7);          % current condition index (0..7), >=8 is condition stop.
        if (icond < 8)
            set(ui.stop_button, 'String', ['Stop Condition ' num2str(icond+1)], 'BackgroundColor', [1 .5 .5]);
        else
            set(ui.stop_button, 'String', 'Condition Stop', 'BackgroundColor', get(ui.save_button, 'BackgroundColor'));
        end
        
        ui.event_summary_index = mod(ui.status_summaries, ui.event_summary_count) + 1;
        ui.status_str = ['nc' num2str(300 + double(ui.status_version(3))) ui.status_version(2) num2str(ui.status_version(1)) '  (Rev ' ui.status_version(4) ')'];
        set(ui.text_version, 'String', ui.status_str);
        set(ui.text_session_info, 'String', ['ID: ' num2str(sid)]);
        
        istream = ui.r(8);        % streaming flag.
        if istream ~= ui.stream_flag
            ui.stream_flag = istream;
            if istream == 0     % Recording stopped
                set(ui.record_button, 'String', 'Record', 'Enable', 'on', 'BackgroundColor', get(ui.save_button, 'BackgroundColor'));
            elseif istream == 1 % Now streaming
                ui.stream_clock = clock;
            end
        end
        if istream == 0
            set(ui.set_button, 'Enable', 'on');
        elseif istream == 1
            set(ui.record_button, 'Enable', 'on', 'BackgroundColor', [0 1 0]);
            set(ui.set_button, 'Enable', 'off');
        elseif istream == 2 % Last start failed
            set(ui.set_button, 'Enable', 'on');
            set(ui.record_button, 'String', 'Record', 'Enable', 'on', 'BackgroundColor', [1 0 0]);
        end
        
        % Grab most recent 12 summaries.  On-going summaries may overwrite
        % anything earlier.
        for summary = ui.status_summaries-11:ui.status_summaries
            index = mod(summary, ui.event_summary_count) + 1;
            offset = double(bitand(summary, 15)); % offset for summary data in ui.r
            imin = (20:2:34) + (offset * 16);
            imax = imin + 256;
            ievt = (533:540) + (offset * p.event_count);
            istm = (661:666) + (offset * p.stim_channel_count);

            for ich=1:p.event_count
                ui.event_summary_min(index,ich) = GetWordReg(imin(ich));
                ui.event_summary_max(index,ich) = GetWordReg(imax(ich));
                ui.event_summary_event(index,ich) = ui.r(ievt(ich));
                if ich <= p.stim_channel_count
                    ui.event_summary_stim(index,ich) = ui.r(istm(ich));
                end
            end
        end;
        
        %Impedance check values
        menu_item = getnum(ui.z_check_menu);
        if menu_item > 0
            menu = get(ui.z_check_menu, 'String');
            z_text = sprintf('Impedance at 1 kHz (%s)\rSubject: %s\rDate: %s\r', menu{menu_item+1},p.subject, datestr(now));
            z_polarity = ' ';
            for ichan = 1:32
                z_kOhm = round(double(GetULongReg(752 + ichan * 4)) / 1000.0);
                if p.intan_channels == 16
                    if z_polarity ~= 'P'  % Differential Intan chips have 16 Pos and 16 Neg input pins.
                        z_polarity = 'P';
                    else
                        z_polarity = 'N';
                    end
                    jchan = floor((ichan + 1) / 2);
                else
                    jchan = ichan;
                end
                z_scale_max = [125, 1250, 12500]; % Maximum of 125 kOhm at the 100 kOhm scale, 1250kOhm at 1MOhm scale, 12500 kOhm at 10MOhm scale.
                if (z_kOhm > z_scale_max(menu_item))
                    z_text = sprintf('%s%s%s (%s) off scale\r', z_text, ui.channel_names{jchan}, z_polarity, p.channel_names{jchan});
                else
                    z_text = sprintf('%s%s%s (%s) %d kOhm\r', z_text, ui.channel_names{jchan}, z_polarity, p.channel_names{jchan}, z_kOhm);
                end
            end
            set(ui.z_check_output, 'String', z_text);
        end
        
        %Test values from ARM chip
        test1 = GetULongReg(1012);
        test2 = GetULongReg(1016);
        if (ui.test1 ~= test1) || (ui.test2 ~= test2)
            disp(['test1 0x' dec2hex(test1) ', test2 0x' dec2hex(test2)]);
            ui.test1 = test1;
            ui.test2 = test2;
        end
    end

    function GetSweepBlock()
        index = GetULongReg(8); % Index of final data point.
        for isample = 1:500
            index = index + 1;  % Index of next data point ...
            if index >= 500      %  ... with wrap around.
                index = 0;
            end
            ui.current_sweep(isample) = GetWordReg(20 + index*2);
        end
        
        ui.last_sweep_index = ui.last_sweep_index + 1;
        if ui.last_sweep_index > ui.max_sweeps
            ui.last_sweep_index = 1;
        end
        ui.sweep_ok(ui.last_sweep_index) = 1;
        ui.sweeps(ui.last_sweep_index, :) = ui.current_sweep;
        %disp(['Sweep ' num2str(ui.last_sweep_index)]);
        %disp([num2str(min(ui.current_sweep)), ', ' num2str(max(ui.current_sweep))]);
    end

    %%%
    % Main communications routine
    function check_send_command()
        % Read last response from neurochip.
        port = irport();
        % disp(['ConStatus ' num2str(port.BytesAvailable)]);
        if port.BytesAvailable < 5
            pause(0.01); % extra time to finish the block.
            if port.BytesAvailable < 5
                ui.irport_status = 3; % Error status on Connect button.
            end
        end

        while port.BytesAvailable >= 5
            % Look for command byte in buffered characters
            try
                byte = fread(port, 1, 'uint8');
            catch 
                byte = [];
            end
            while ~isempty(byte)
                if port.BytesAvailable < 4
                    break;
                end
                if byte == 'A'
                    % Response packet from FPGA SDCard Buffer Address Read.
                    try
                        %disp('A');
                        val = fread(port, 4, 'uint8');
                        fpga_int = bitor(bitshift(val(1),24), bitor(bitshift(val(2),16), bitor(bitshift(val(3), 8), val(4))));
                        set(ui.text_record_info, 'String', ['SD: ' num2str(fpga_int)]);
                    catch
                    end
                    byte = [];
                elseif ui.command_get && (byte == 'p')
                    % Get Parameter block 4 byte header with 2044 byte following data.
                    ui.r(1) = byte;
                    for ibyte = 2:4  % Requires 3 packet_marker bytes in a row.
                        try
                            byte = fread(port, 1, 'uint8');
                        catch
                            byte = [];
                        end
                        if byte ~= ui.packet_marker
                            break;
                        end
                        ui.r(ibyte) = byte;
                    end
                    for i = 1:20  % Wait for entire packet to arrive
                        if port.BytesAvailable >= ui.registerBlockSize-4
                            break;
                        end
                        pause(0.05);
                    end
                    if (byte ~= ui.packet_marker) || (port.BytesAvailable < ui.registerBlockSize-4)
                        ui.irport_status = 2; % Bad packet header warning.
                        if (byte == ui.packet_marker)
                            % Correct header, assume lost packet.
                            byte = [];
                            try
                                fread(port, port.BytesAvailable, 'uint8');
                            catch
                            end
                        end
                    else
                        % Read data block
                        try
                            ui.r(5:ui.registerBlockSize) = fread(port, ui.registerBlockSize-4, 'uint8');
                            checksum = sum(ui.r(5:ui.registerBlockSize-4));
                        catch
                            checksum = bitcmp(GetULongReg(ui.registerBlockSize-4)); % make a bad checksum
                        end
                        if checksum ~= GetULongReg(ui.registerBlockSize-4)
                            ui.irport_status = 2;  % Bad checksum warning.
                        else
                            GetParameterBlock();
                            ui.command_get = 0;    % Get command finished OK.
                            ui.irport_status = 1;  % Clear IR status.
                        end
                        byte = []; % Done with this block
                    end
                elseif (byte == 'i') || (byte == 't')
                    % Status block 4 byte header with 1020 byte following data.
                    ui.r(1) = byte;
                    %disp(['Command ' byte]);
                    for ibyte = 2:4  % Requires 3 packet_marker bytes in a row.
                        try
                            byte = fread(port, 1, 'uint8');
                        catch
                            byte = [];
                        end
                        if byte ~= ui.packet_marker
                            break;
                        end
                        ui.r(ibyte) = byte;
                    end
                    for i = 1:10  % Wait for entire packet to arrive
                        if port.BytesAvailable >= 1020
                            break;
                        end
                        pause(0.05);
                    end
                    if (byte ~= ui.packet_marker) || (port.BytesAvailable < 1020)
                        ui.irport_status = 2; % Bad packet warning.
                        if (byte == ui.packet_marker)
                            % Correct header, assume lost packet.
                            byte = [];
                            try
                                fread(port, port.BytesAvailable, 'uint8');
                            catch
                            end
                        end
                    else
                        % Read data block
                        try
                            %disp(['before fread ' num2str(port.BytesAvailable)]);
                            ui.r(5:1024) = fread(port, 1020, 'uint8');
                            %disp(['after fread ' num2str(port.BytesAvailable)]);
                        catch
                            SetULongReg(1020,0); % Bad data block.
                        end
                        if ui.packet_check ~= GetULongReg(1020)
                            ui.irport_status = 2;  % Bad checksum warning.
                            %disp('Bad Packet Check');
                        else
                            if (ui.r(1) == 'i')
                                GetStatusBlock();
                                if ui.parameter_set_status == 1
                                    ui.command_set = 0; % Set command finished OK.
                                    HighlightSetButton(0);
                                end
                                ui.irport_status = 1;  % Clear IR status.
                            else
                                GetSweepBlock(); % This is a sweep block.
                            end  
                        end
                        byte = []; % Done with this block
                    end
                else
                    byte = [];
                end
            end % While ~isempy(byte)
        end % While port

        % Clear any bytes left over in the read buffer.
        if port.BytesAvailable > 0
            ui.irport_status = 2; % Warning, discarded bytes.
            try
                fread(port, port.BytesAvailable, 'uint8');
            catch
            end
        end
        
        % Check for pending commands
        if (ui.command == 's') || ui.command_set
            % Send current parameters.  NC3 will respond with a status block.
            %disp(['ui.command ' ui.command]);
            ui.command = 0;
            ui.command_set = 1;    % Set parameters command is pending.
            ui.command_set = 0;    % ^^^TODO remove this test
            ui.irport_status = 2;  % Waiting for parameters.
            SetParameterBlock();
            % Apparently versions of Matlab can't write much more than 1k
            % bytes at a time to a serial port, so we break the write up.
            for startbyte = 1:1024:ui.registerBlockSize
                endbyte = min(startbyte+1023, ui.registerBlockSize);
                while strcmp(port.TransferStatus, 'write')
                    pause(0.01);
                end
                fwrite(port, ui.r(startbyte:endbyte), 'uint8', 'async');
            end              
        elseif (ui.command == 'g') || ui.command_get
            % Request current parameters.  NC3 will respond with a
            % parameter block and a status block.
            ui.command = 0;
            ui.command_get = 1;   % Get parameters command is pending.
            ui.irport_status = 2; % Waiting for request.
            ir_write(uint8(['g', ui.packet_marker, ui.packet_marker, ui.packet_marker]));
        else
            % Send status info command. NC3 will respond with a status block and
            % an optional sweep block, and a read of the SDCard Current Block Address.
            filter_chan = get(ui.event_list, 'Value') - 1;  % Selected filter channel 0..7
            %ir_write(uint8(['i', sweep_decimation(), filter_chan, ui.packet_marker]));
            ir_write(uint8('R'), uint16(24582), uint32(0), uint8(['i', sweep_decimation(), filter_chan, ui.packet_marker]));
            %disp(['txstatus ' num2str(toc)]);
            ui.command = 0;            
        end
        
        % Update the record button with elapsed streaming time information.
        if ui.stream_flag == 1
            elapsed_time = etime(clock, ui.stream_clock);
            hours = floor(elapsed_time / 3600);
            elapsed_time = elapsed_time - (hours * 3600);
            mins = floor(elapsed_time / 60);
            elapsed_str = sprintf('%02d:%02d:%02d', hours, mins, floor(elapsed_time - (mins * 60)));
            set(ui.record_button, 'String', elapsed_str);
        end
    end

    function decimation = sweep_decimation()
        % Calculate decimation for sweep data from the neurochip.
        % Sweeps are currently 500 samples with 100 pre-event, and 400
        % post-event samples.
        curr_event_chan = get(ui.event_list, 'Value'); % Selected Event channel
        if ismember(p.event_source1(curr_event_chan), ui.channel_Accel_index) || ismember(p.event_source2(curr_event_chan), ui.channel_Accel_index)
            decimation = 1; % Accelerometer data doesn't decimate.
        else
            xsamples = ui.ms2fs * getnum(ui.ms_scale);
            if xsamples <= 0
                xsamples = p.event_window2_delay(curr_event_chan);
                if xsamples <= 0
                    xsamples = 400;
                end
            end
            decimation = max(1, min(255, floor((xsamples-1) / 400) + 1));
        end
    end

    function update_all_controls()
        % match control values to current parameter values.
        RecordSelect(0,0);
        ConditionSelect(0,0);
        EventSelect(0,0);
        for i=p.stim_req_count:-1:1
            set(ui.req_list, 'Value', i);
            ReqSelect(0,0);
        end
        TestSelect(0,0);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nested Callback functions


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Handle a click in the Menu Bar Button
    function MenuBarButton(object, event) %#ok<INUSD>
        ui.sweep_ok(:) = 0;
        on_off = get(ui.figure, 'MenuBar');
        switch on_off
            case 'none'
                set(ui.figure,'MenuBar','figure','Toolbar','Figure');
            case 'figure'
                set(ui.figure,'MenuBar','None','Toolbar','None');
        end % end switch/case
    end % end MenuBarButton nested function


    function SaveFolderButton(object, event) %#ok<INUSD>
        if ui.save_in_progress
            return;
        end
        ct = datevec(now);
        fileindex = 1;
        subject = get(ui.subject, 'string');
        datestamp = [num2str(ct(1)) num2str(ct(2),'%02d') num2str(ct(3),'%02d') '_' num2str(fileindex,'%02d')];
        foldername = sprintf('%s\\%s\\%s_%s', ui.datapath, subject, subject, datestamp);
        while exist(foldername, 'file')
            fileindex = fileindex + 1;
            datestamp = [num2str(ct(1)) num2str(ct(2),'%02d') num2str(ct(3),'%02d') '_' num2str(fileindex,'%02d')];
            foldername = sprintf('%s\\%s\\%s_%s', ui.datapath, subject, subject, datestamp);
        end
        
        options.Resize = 'on';
        options.WindowStyle = 'modal';
        options.Interpreter = 'none';
        a = inputdlg({'SDCard drive:', 'Save in folder:'},  'Download SDCard data to folder', [1 100], {ui.sddrive, foldername}, options);

        if ~isempty(a)
            foldername = a{2};
            ui.sddrive = [upper(deblank(a{1})) 'E'];
            ui.sddrive = ui.sddrive(1);
            if (ui.sddrive < 'D') || (ui.sddrive > 'Z')
                ui.sddrive = 'E';
            end
            
            % Create a data folder to hold all of the supporting files.
            while (foldername(end) == '\') || (foldername(end) == '/')
                foldername = foldername(1:end-1);
            end
            if ~exist(foldername, 'dir')
                [status,message,messageid] = mkdir(foldername); %#ok<ASGLU>
                if status == 0
                    warndlg('Download could not create the data folder','Save Error');
                    return;
                end
            else
                warndlg('Download cannot overwrite an existing data folder','Save Error');
                return;
            end 
            
            % File prefix is a copy of the last folder name.
            index = find((foldername == '\') | (foldername == '/'));
            if isempty(index)
                destprefix = [foldername '\' foldername];
            else
                destprefix = [foldername foldername(index(end):end)];
            end
           
            
            ui.save_in_progress = 1;
            set(ui.savefolder_button, 'BackgroundColor', [1 1 0]);
            rate(1:35) = p.channel_rate(1:35); % Limit to Intan channels
            drawnow;
            id = nc323save(ui.sddrive, 0, 172800, destprefix, rate, ui.savefolder_button, 'String', '00:00:00');
            set(ui.savefolder_button, 'BackgroundColor', get(ui.save_button, 'BackgroundColor'));
            ui.save_in_progress = 0;
            
            if (id(1) < 0)
                disp(['Error ' num2str(-id(1)) ': Session ID not found.  File on SDCard may be corrupt.']);
                warndlg('Could not find Session ID, file on SDCArd may be corrupt.','Save Error');
            else
                % Save settings to a file.
                subject = deblank(get(ui.subject, 'String'));
                subject(subject == ' ') = '_';       
                settings_folder = [ui.path subject];
                if ~exist(settings_folder, 'dir')
                    mkdir(settings_folder);
                end          
                file_name = sprintf('%s\\session_%.5d.mat', settings_folder, id(1)); 
                try
                    copyfile(file_name, [destprefix '.mat'], 'f');
                catch
                    warndlg(['Warning: Please copy settings for session ' num2str(id(1)) ' to ' destprefix '.mat']);
                end
            end
            
            % Save a text info file describing contents.
            try
                fid = fopen([destprefix '_info.txt'], 'w');
                if fid ~= -1
                    fprintf(fid, 'Neurochip3 Data\r\n');
                    fprintf(fid, 'Version %d\r\n', ui.version);
                    fprintf(fid, 'Session %d\r\n', id(1));
                    fprintf(fid, 'Name %s\r\n', foldername);
                    fprintf(fid, 'Date %s\r\n', date);
                    fprintf(fid, 'File _ChanN.i16 contains int16 samples for analog Chan N Rate Suffix\r\n');
                    for ichan=1:35
                        if id(ichan+1) > 0
                            if (ichan < 33)
                                fprintf(fid, 'Chan %d %d _Chan%02d.i16\r\n', ichan, id(ichan+1), ichan);
                            else
                                fprintf(fid, 'Chan %d %d _Chan%02d.u16\r\n', ichan, id(ichan+1), ichan);
                            end
                        end
                    end
                    fprintf(fid, 'Chan 36 100 _AccelM.i16\r\n');
                    fprintf(fid, 'Chan 37 100 _AccelX.i16\r\n');
                    fprintf(fid, 'Chan 38 100 _AccelY.i16\r\n');
                    fprintf(fid, 'Chan 39 100 _AccelZ.i16\r\n');
                    fprintf(fid, 'Chan 40 100 _AccelT.i16\r\n');
                    fprintf(fid, 'File _Events.u32 contains (ID, Value, Timestamp5k) uint32 triplets for each event\r\n');
                    fprintf(fid, 'Event ID 0 = Condition_Start, Value = condition index 1..9\r\n');
                    fprintf(fid, 'Event ID 1 = Discrimination, Value = discriminator index 1..8\r\n');
                    fprintf(fid, 'Event ID 2 = Stimulus, Value = stimulus index 1..6\r\n');
                    fprintf(fid, 'Event ID 3 = Sampling_Status, Value = 0 (Samples stopped), 1 (Samples started)\r\n');
                    fprintf(fid, 'Event ID 4 = One_Second_Marker, Value = packet index\r\n');
                    %fprintf(fid, 'Event ID 5 = CRC_Error, Value = (Calculated_CRC << 16) | Expected_CRC\r\n');
                    fprintf(fid, 'File _Digi00.u16 contains uint16 header word for each sample packet (5000 packets per second)\r\n');
                    fprintf(fid, 'Diginfo hexadecimal bit fields in _Digi00.u16\r\n');                    
                    fprintf(fid, 'Diginfo 0x8000 Bad Packet Check\r\n');
                    fprintf(fid, 'Diginfo 0x4000 One Second Marker\r\n');
                    fprintf(fid, 'Diginfo 0x2000 Analog Samples Included\r\n');
                    fprintf(fid, 'Diginfo 0x1000 Stimulus Trigger Bit\r\n');
                    fprintf(fid, 'Diginfo 0x0F00 Stimulus Channel Index = 1+bitand(bitshift(diginfo, -8), 15);\r\n');
                    fprintf(fid, 'Diginfo 0x0080..0x0001 Discriminator 8..1 Trigger Bits\r\n');
                    fclose(fid);
                end
            catch
                fclose(fid);
                disp(['Warning: Text information file is incomplete for ' foldername]);
            end                
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Handle a click in the Test Button
    %This button is normally hidden and is only used for testing purposes.
    function TestButton(object, event) %#ok<INUSD>
    % Currently this button dumps the stimulus sequence tables to the
    % console.
    
        % Pause the Connect timer so we can squeeze in a test command.   
        if (ui.connection_status)
            stop(ui.command_timer);
            pause(0.5);
        end

        % For reads, clear any bytes left over in the read buffer.
        port = irport();
        if (port.BytesAvailable > 0)
            fread(port, port.BytesAvailable, 'uint8');
        end

        % Read the ADC table.  
        disp('ADC TABLE START 0x0200');
        addr = hex2dec('0200');
        ch = 'R';

        for icount = 1:64
            ir_write(uint8(ch), uint16(addr), uint32(0));
          
            % Wait for reponse from a read request.
            
            pause(.25);
            if port.BytesAvailable >= 5
                % Look for command byte in buffered characters
                byte = fread(port, 1, 'uint8');
                if byte == 'A'
                    val = fread(port, 4, 'uint8');
                    fpga_int = bitor(bitshift(val(1),24), bitor(bitshift(val(2),16), bitor(bitshift(val(3), 8), val(4))));
                    if fpga_int < 65536
                        disp(['0x' dec2hex(addr, 4) ' 0x' dec2hex(fpga_int, 4)]);
                    else
                        disp(['0x' dec2hex(addr, 4) ' 0x' dec2hex(fpga_int, 8)]);
                    end
                end
            end

            % Clear any bytes left over in the read buffer.
            if port.BytesAvailable > 0
                fread(port, port.BytesAvailable, 'uint8');
            end

            addr = addr + 1;
        end
        
%         % Read stim table.
%         disp('SEQUENCER START 0x1400');
%         addr = hex2dec('1400');
%         ch = 'R';
% 
%         for icount = 1:128
%             ir_write(uint8(ch), uint16(addr), uint32(0));
%           
%             % Wait for reponse from a read request.
%             
%             pause(.25);
%             if port.BytesAvailable >= 5
%                 % Look for command byte in buffered characters
%                 byte = fread(port, 1, 'uint8');
%                 if byte == 'A'
%                     val = fread(port, 4, 'uint8');
%                     fpga_int = bitor(bitshift(val(1),24), bitor(bitshift(val(2),16), bitor(bitshift(val(3), 8), val(4))));
%                     if fpga_int < 65536
%                         disp(['0x' dec2hex(addr, 4) ' 0x' dec2hex(fpga_int, 4)]);
%                     else
%                         disp(['0x' dec2hex(addr, 4) ' 0x' dec2hex(fpga_int, 8)]);
%                     end
%                 end
%             end
% 
%             % Clear any bytes left over in the read buffer.
%             if port.BytesAvailable > 0
%                 fread(port, port.BytesAvailable, 'uint8');
%             end
% 
%             addr = addr + 1;
%         end
%   
%         % Read stim table.  
%         disp('SEQUENCER TABLE 0x1800');
%         addr = hex2dec('1800');
%         ch = 'R';
% 
%         for icount = 1:256
%             ir_write(uint8(ch), uint16(addr), uint32(0));
%           
%             % Wait for reponse from a read request.
%             
%             pause(.25);
%             if port.BytesAvailable >= 5
%                 % Look for command byte in buffered characters
%                 byte = fread(port, 1, 'uint8');
%                 if byte == 'A'
%                     val = fread(port, 4, 'uint8');
%                     fpga_int = bitor(bitshift(val(1),24), bitor(bitshift(val(2),16), bitor(bitshift(val(3), 8), val(4))));
%                     if fpga_int < 65536
%                         disp(['0x' dec2hex(addr, 4) ' 0x' dec2hex(fpga_int, 4)]);
%                     else
%                         disp(['0x' dec2hex(addr, 4) ' 0x' dec2hex(fpga_int, 8)]);
%                     end
%                 end
%             end
% 
%             % Clear any bytes left over in the read buffer.
%             if port.BytesAvailable > 0
%                 fread(port, port.BytesAvailable, 'uint8');
%             end
% 
%             addr = addr + 1;
%         end

        % Restart the Connect timer if necessary
        if (ui.connection_status)
            start(ui.command_timer);  
        end        
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Handle a click in the Clear Sweep Button
    function ClearButton(object, event) %#ok<INUSD>
        ui.sweep_ok(:) = 0;
        update_plots();
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Handle a click in the Clear Stop Condition Button
    %This button sends an IR packet to the Neurochip3 in an attempt
    %to stop all stimulation by entering the "Stop Condition".
    function StopButton(object, event) %#ok<INUSD>
        stop(ui.command_timer);  % Stop connect timer.
        set(ui.stop_button, 'String', 'Stopping...');
        pause(1);                % Wait for possible packet in progress.
        ir_write(uint8(['h', 165, 165, 165]));
        pause(0.5);              % Do it again to be sure
        ir_write(uint8(['h', 165, 165, 165]));
        pause(0.5);              % Do it again to be sure
        ir_write(uint8(['h', 165, 165, 165]));
        start(ui.command_timer);  % Start connect timer.
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Handle a click in the Connect Button.  This starts the command_timer
    %which periodically calls the CommandTimerFcn.
    
    function CommandTimerFcn(object, event) %#ok<INUSD>
        % Check for received data and send requested commands to the neurochip
        check_send_command();
        % Handle animations
        set(ui.connect_button, 'BackgroundColor', ui.status_color{ui.irport_status});
        update_plots();
    end

    function ConnectButton(object, event) %#ok<INUSD>
        if strcmp(get(ui.connect_button, 'String'), 'Connect')
            irport_init();
            if ~isempty(irport())
                set(ui.connect_button, 'String', 'Disconnect');
                set(ui.record_button, 'Enable', 'on', 'BackgroundColor', get(ui.save_button, 'BackgroundColor'));
                set(ui.comport_popup, 'Enable', 'off');
                set(ui.set_button, 'Enable', 'on');
                set(ui.stop_button, 'Enable', 'on');
                ui.connection_status = 1;
                ui.irport_status = 2;
                ui.command = 0;      % Clear any pending command.
                ui.command_set = 0;
                ui.command_get = 0;
                port = irport();
                if port.BytesAvailable > 0
                    try % Clear any bytes in serial port buffer.
                        fread(port, port.BytesAvailable, 'uint8');
                    catch
                    end
                end
                start(ui.command_timer);
            else
                beep;
            end
        else
            ui.connection_status = 0;
            stop(ui.command_timer);
            set(ui.record_button, 'Enable', 'off', 'String', 'Record');
            set(ui.comport_popup, 'Enable', 'on');
            set(ui.open_button, 'Enable', 'on');
            set(ui.connect_button, 'String', 'Connect', 'BackgroundColor', get(ui.save_button, 'BackgroundColor'));
            set(ui.set_button, 'Enable', 'off');
            set(ui.stop_button, 'Enable', 'off');
        end
    end % end ConnectButton nested function


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Handle files 
    
    function load_settings(file_name)
        % Load settings from saved parameter files.
        try
            par = load(file_name, 'p');
            if isfield(par, 'p')
                par = par.p;
                names = fieldnames(par);
                n = length(names);
                for i=1:n
                    p.(names{i}) = par.(names{i});
                end
                if p.version < ui.version
                    if p.version < 310
                        % Convert old channel name list to new list.
                        p.intan_channels = 16;
                        p.remote_channels = 0;
                        p.channel_count = 43;
                        p.channel_names{33} = p.channel_names{17}; % Move Aux channels higher in list
                        p.channel_names{34} = p.channel_names{18};
                        p.channel_names{35} = p.channel_names{19};
                        for i=17:ui.channel_count % Extra channels for 32 channel Intan chip
                            p.channel_names{i} = ui.channel_names{i};
                        end
                        
                        % Fix up channel rates to match new name order
                        p.channel_rate(33:35) = p.channel_rate(17:19); 
                        p.channel_rate(17:32) = 0;    % Old files only had first 16 channels
                        p.channel_rate(36:40) = 100;  % Accelerometer channels
                        p.channel_rate(41:48) = 2000; % Remote emg channels
                    end
                    
                    if p.version < 320
                        % Version 320 increase conditions to 8 and stimulus channels to 6.
                        % Fill in missing values with defaults.
                        p.stim_cond_count = 8; % Number of possible stimulus conditions is now 8
                        p.stim_channel_count = 6; % Number of stimulus channels is now 6
                        p.event_refractory_period = repmat(5, p.event_count, 1);   % Refractory period when refractory_count reached.
                        p.event_refractory_count = ones(p.event_count, 1);    % Count triggering refractory period.
                        p.event_refractory_window = repmat(5, p.event_count, 1);   % Counts must all happen inside of given millisecond window.
                        p.stim_cond_duration(5:8) = 0;  % Condition duration in seconds, 0 for unlimited duration.
                        p.stim_cond_repeat_condition(5:8) = [5 6 7 8]; % Go to this condition (zero based index) when duration expires.
                        p.stim_cond_repetitions(5:8) = 0;  % Maximum number of repetions before continuing to following condition
                        p.stim_cond_stim_mode(5:8) = 0;    % 0=No stimulation, Otherwise this is the number of active stim channels.
                        for icond = 1:8
                            for istim = 1:6
                                if (icond >= 5) || (istim >= 4)
                                    p.stim_isi(icond, istim) = 2;      % Minimum Inter-Stimulus Interval
                                    p.stim_level(icond, istim) = 0;    % Stimulation level in uA.
                                    p.stim_width(icond, istim) = 0.2;  % Phase width in ms.
                                    p.stim_delay(icond, istim) = 0;    % Millisecond Stimulus delay from trigger.
                                    p.stim_pulses(icond, istim) = 0;   % Pulses in stimulus pulse train. 0 = Don't deliver or mark the stimulus event.
                                    p.stim_level2(icond, istim) = 0;   % Optional level and width for the second phase (of the bi-phasic pulse)
                                    p.stim_width2(icond, istim) = 0.2;
                                end 
                                % Defaults for parameters added in nc320.
                                p.stim_cathodic_first(icond, istim) = 0; % Pin# for Positive voltage on first phase.  Positive current flows out of this pin on the first phase
                                p.stim_anodic_first(icond, istim) = 0;   % Pin# for Negative voltage on first phase. 0 
                                p.stim_limit_count(icond, istim) = 200;  % Safe stimulation limit. Exceeding this causes a time-out.
                                p.stim_limit_window(icond, istim) = 1;   % Safe as long as stimulations <= count inside each window (given in seconds)
                                p.stim_limit_timeout(icond, istim) = 5;  % Number of seconds to stop stimulation if safety limit exceeded.
                            end
                        end
                    end
                    
                    % Parameters to fix up for all versions.
                    p.stim_fast_settle_menu = get(ui.stim_fast_settle, 'String');
                    p.channel_rate(36:40) = 100;  % Accelerometer channels
                    p.channel_rate(41:48) = 2000; % Remote emg channels
                    p.version = ui.version; % Update version number to current
                end
                for i=36:40
                    p.channel_names{i} = ui.channel_names{i}; % Accelerometer channel names cannot be changed.
                end
                set(ui.settings_file, 'String', file_name);
            end
            update_all_controls();
        catch
            beep;
        end
    end

    function load_last_settings()
        % Called during initialzation to load the last saved settings file
        if exist(ui.last_settings, 'file')
            fid = fopen(ui.last_settings, 'r');
            if fid ~= -1
                file_name = fgetl(fid);
                comportstr = fgetl(fid);
                comportnum = min(12, max(1, str2double(comportstr(4:end))));
                ui.session_id = str2double(fgetl(fid));
                ui.sddrive = [fgetl(fid) 'E'];
                ui.sddrive = ui.sddrive(1);
                set(ui.comport_popup, 'Value', comportnum);
                fclose(fid);
                load_settings(file_name)
            end
        end
    end

    function save_settings(file_name)
        try
            fid = -1;
            p.date = date;
            p.neurochip_firmware = get(ui.text_version, 'String');
            save(file_name, 'p');
            if nargin == 1
                set(ui.settings_file, 'String', file_name);
                fid = fopen(ui.last_settings, 'w');
                comport = min(10, max(1, get(ui.comport_popup, 'Value')));
                fprintf(fid, '%s\r\nCom%d\r\n%d\r\n%s\r\n', file_name, comport, ui.session_id, ui.sddrive);
                fclose(fid);
            end
        catch
            if fid ~= -1
                fclose(fid);
            end
        end
    end

    %%%
    %Handle a click in the Record Button
    function RecordButton(object, event) %#ok<INUSD>
        port = irport();
        set(ui.record_button, 'BackgroundColor', [1 1 0]);
        set(ui.set_button, 'Enable', 'off');
        stop(ui.command_timer);  % Stop connect timer.
        pause(1);                % Wait for possible packet in progress.
        
        if ~strcmp(get(ui.record_button, 'String'), 'Record')
            % Stop recording, close the SDCard.
            ir_write(uint8(['c', 165, 165, 165]));
        else
            % Start recording
            ir_write(uint8(['o', 165, 165, 165]));

            % Save current settings to a file.
            p.session_id = ui.session_id;
            p.subject = deblank(get(ui.subject, 'String'));
            p.subject(p.subject == ' ') = '_';       
            settings_folder = [ui.path p.subject];
            if ~exist(settings_folder, 'dir')
                mkdir(settings_folder);
            end          
            file_name = sprintf('%s\\session_%.5d.mat', settings_folder, p.session_id); 
            save_settings(file_name);         
        end
        
        % Wait a bit for command to finish then clear input buffer and restart Connect updates
        pause(2);
        if port.BytesAvailable > 0
            try % Clear any bytes in serial port buffer.
                fread(port, port.BytesAvailable, 'uint8');
            catch
            end
        end
        start(ui.command_timer);
    end % end RecordButton nested function

    %%%
    %Handle a click in the Set Button
    function SetButton(object, event) %#ok<INUSD>
        ui.command = 's';
        ui.session_id = ui.session_id + 1; % increment session ID.
        if (ui.session_id >= 65536)
            ui.session_id = 1;
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Handle a click in the Open Button
    function OpenButton(object, event) %#ok<INUSD>
        [filename, pathname] = uigetfile('*.nc3;*.mat', 'Select Neurochip settings file (.nc3) file');
        if isequal(filename,0) || isequal(pathname,0)
            % disp('User pressed cancel')
        else
            load_settings([pathname filename]);
            ui.command = 's'; % Update neurochip if connected.
        end
    end % end OpenButton nested function

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Handle a click in the Save Button
    function SaveButton(object, event) %#ok<INUSD>
        [filename, pathname] = uiputfile('*.nc3;*.mat', 'Save Neurochip settings file (.nc3) file');
        if isequal(filename,0) || isequal(pathname,0)
            % disp('User pressed cancel')
        else
            p.filename = [pathname filename];
            save_settings(p.filename);
        end
    end % end OpenButton nested function

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Handle a click in one of the tab buttons
    function TabSelect(object, event) %#ok<INUSD>
        for i = 1:length(ui.tab_names)
            if object == ui.tab_buttons(i)
                panel_color = get(ui.tab_panels(i), 'BackgroundColor');
                pos = get(ui.tab_buttons(i), 'Position');
                pos = [pos(1) + 0.02, pos(2) - 0.05, pos(3) - 0.06, 0.1];
                set(ui.tab_buttons(i), 'BackgroundColor', panel_color);
                set(ui.tab_cover, 'BackgroundColor', panel_color, ...
                    'ForegroundColor', panel_color, 'Position', pos);
                set(ui.tab_panels(i), 'Visible', 'on');
            else
                set(ui.tab_buttons(i), 'BackgroundColor', get(ui.figure, 'Color'));
                set(ui.tab_panels(i), 'Visible', 'Off');
            end
            if getnum(ui.z_check_menu) > 0
                % Warn user that impedance checking is on with yellow tab name.
                set(ui.tab_buttons(4), 'BackgroundColor', 'y');
            end
        end
    end

    %%%
    % User clicked in one of the test buttons, send test over the Irport.
    function TestFPGAButton(object, event) %#ok<INUSD>
        index = find(ui.test_fpga_button == object);
        if ~isempty(index) && (index <= 2*p.test_count)
            port = irport();
            ui.last_test_index = index;
            if index > p.test_count
                ch = 'R'; % Read button pressed
                index = index - p.test_count;
                set(ui.test_response(index), 'String', ''); % Blank out last response.
            else
                ch = 'W'; % Write button pressed
            end
            
            % Pause the Connect timer so we can squeeze in a test command.   
            if (ui.connection_status)
                stop(ui.command_timer);
                pause(0.5);
            end
            
            % For reads, clear any bytes left over in the read buffer.
            if (ch == 'R') && (port.BytesAvailable > 0)
                fread(port, port.BytesAvailable, 'uint8');
            end

            % Send the command                
            ir_write(uint8(ch), uint16(p.test_address(index)), uint32(p.test_value(index)));
          
            % Wait for reponse from a read request.
            if (ch == 'R')
                pause(.25);
                if port.BytesAvailable >= 5
                    % Look for command byte in buffered characters
                    byte = fread(port, 1, 'uint8');
                    if byte == 'A'
                        val = fread(port, 4, 'uint8');
                        fpga_int = bitor(bitshift(val(1),24), bitor(bitshift(val(2),16), bitor(bitshift(val(3), 8), val(4))));
                        set(ui.test_response(index), 'String', ['0x' dec2hex(fpga_int, 8)]);
                    end
                else
                    %disp(num2str(port.BytesAvailable));
                end
                % Clear any bytes left over in the read buffer.
                if port.BytesAvailable > 0
                    fread(port, port.BytesAvailable, 'uint8');
                end
            end
            
            % Restart the Connect timer if necessary
            if (ui.connection_status)
                start(ui.command_timer);  
            end        
        end
    end

    function TestParam(object, event) %#ok<INUSD>
        % UI test parameter changed.  Convert them to p.test parameters
        for i=1:p.test_count
            p.test_name{i} = get(ui.test_name(i), 'String');
            p.test_address(i) = getnum(ui.test_address(i));
            p.test_value(i) = getnum(ui.test_value(i));
        end
    end

    function TestSelect(object, event) %#ok<INUSD>
        % Copy all p.test parameters to the user interface.
        for i=1:p.test_count
            set(ui.test_name(i), 'String',  p.test_name{i});
            set(ui.test_address(i), 'String', ['0x' dec2hex(p.test_address(i), 4)]);
            set(ui.test_value(i), 'String', ['0x' dec2hex(p.test_value(i), 8)]);
        end
    end

    %%%
    % Update GUI to reflect changes in recording parameters
    function RecordSelect(object, event) %#ok<INUSD>
        setnum(ui.upper_bandwidth_popup, p.analog_upper_bandwidth_code);
        setnum(ui.lower_bandwidth_popup, p.analog_lower_bandwidth_code);
        setnum(ui.dsp_high_pass_popup, p.dsp_high_pass_code);
        setnum(ui.intan_channels_popup, (p.intan_channels / 16) - 1);
        setnum(ui.remote_channels_popup, p.remote_channels);
        RecordParam(0,0);
        ChanListChanged(0,0);
    end

    %%%
    % Bold/unbold Set Button to indicate parameter changes need to be
    % downloaded to the Neurochip3.
    function HighlightSetButton(flag)
        if (flag == 0)
            set(ui.set_button, 'FontSize', get(ui.save_button, 'FontSize'), 'FontWeight', 'normal');
        else
            set(ui.set_button, 'FontSize', 12, 'FontWeight', 'bold');
        end
    end

    %%%
    %Handle a changed recording parameter
    function RecordParam(object, event) %#ok<INUSD>
        HighlightSetButton(1);
        list = get(ui.upper_bandwidth_popup, 'String');
        p.analog_upper_bandwidth_code = getnum(ui.upper_bandwidth_popup);
        p.analog_upper_bandwidth = str2Hz(list{p.analog_upper_bandwidth_code + 1});
        list = get(ui.lower_bandwidth_popup, 'String');
        p.analog_lower_bandwidth_code = getnum(ui.lower_bandwidth_popup);
        p.analog_lower_bandwidth = str2Hz(list{p.analog_lower_bandwidth_code + 1});
        p.dsp_high_pass_code = getnum(ui.dsp_high_pass_popup);
        p.dsp_high_pass_info = ui.dsp_high_pass_info{p.dsp_high_pass_code + 1};
        p.intan_channels = 16 + 16 * getnum(ui.intan_channels_popup);
        p.remote_channels = getnum(ui.remote_channels_popup);
        
        selected = get(ui.chan_list, 'Value');
        selnames = get(ui.chan_list, 'String');
        n = length(selected);
        name = get(ui.chan_name, 'String');
        rate = ui.chan_rate_conversion(getnum(ui.chan_rate) + 1);
        for i=1:n
            if object == ui.chan_name
                if n == 1
                    p.channel_names{selected(i)} = name;
                else
                    p.channel_names{selected(i)} = [name num2str(i)];
                end
            elseif object == ui.chan_rate
                sel_name = selnames{selected(i)};
                sel_num = str2double(sel_name(1:strfind(sel_name, ')')-1));
                p.channel_rate(sel_num) = rate;
                if ismember(sel_num, ui.channel_Rem_index)
                	p.channel_rate(sel_num) = 2000; % Remote emg channels
                elseif ismember(sel_num, ui.channel_Accel_index)
                	p.channel_rate(sel_num) = 100; % Accelerometer channels
                end
            end
        end
        
        list = {};
        ilist = 0; % index into list
        for i=1:ui.channel_count
            if (i <= p.intan_channels) || ismember(i, ui.channel_Aux_index) || (ismember(i, ui.channel_Rem_index) && (i < ui.channel_Rem_index(1) + p.remote_channels))
                ilist = ilist + 1;
                ratestr = ', Off';
                if p.channel_rate(i) > 0
                    ratestr = [', ' num2str(floor(p.channel_rate(i) / 1000)) ' kHz'];
                end
                try
                    list{ilist} = [num2str(i) ') ' p.channel_names{i} ratestr]; %#ok<AGROW>
                catch
                    p.channel_names{i} = ui.channel_names{i};
                    list{ilist} = [num2str(i) ') ' p.channel_names{i} ratestr]; %#ok<AGROW>
                end
            elseif ismember(i, ui.channel_Accel_index)
                p.channel_rate(i) = 100;
            elseif ismember(i, ui.channel_Rem_index)
                p.channel_rate(i) = 2000;
            else
                p.channel_rate(i) = 0;  % If channel is not in selection list, set its rate to 0.
            end
        end
        selected(selected > ilist) = [];  % Clear invalid selections
        set(ui.chan_list, 'String', list, 'Value', selected);
 
        % Repopulate the event source menues
        
        source_names = {'0'};
        for i=1:16
            % 16 Intan channels
            source_names{end+1} = p.channel_names{i}; %#ok<AGROW>
        end
        if p.intan_channels > 16
            % Channels for the 32 input Intan chip
            for i=17:32
                source_names{end+1} = p.channel_names{i}; %#ok<AGROW>
            end
        end
        for i=33:40
            % Aux and Accel channels
            source_names{end+1} = p.channel_names{i}; %#ok<AGROW>
        end
        if p.remote_channels > 0
            % Remote emg channels
            for i=41:40+p.remote_channels
                source_names{end+1} = p.channel_names{i}; %#ok<AGROW>
            end
        end

        set(ui.event_source1, 'String', source_names);
        selected = get(ui.event_source1, 'Value');
        if selected > length(source_names)
            set(ui.event_source1, 'Value', 1);
        end
        
        set(ui.event_source2, 'String', source_names);
        selected = get(ui.event_source2, 'Value');
        if selected > length(source_names)
            set(ui.event_source2, 'Value', 1);
        end
         
        % Check for changes in sample rates that could affect filter coeffs.
        for event_index = 1:8
            calculate_filter_coefficients(event_index);
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Handle a change in the channel selection list
    function ChanListChanged(object, event) %#ok<INUSD>
        selected = get(ui.chan_list, 'Value');
        if isempty(selected)
            set(ui.chan_name, 'Enable', 'Off');
            set(ui.chan_rate, 'Enable', 'Off');
        else            
            selnames = get(ui.chan_list, 'String');
            sel_name = selnames{selected(1)};
            sel_num = str2double(sel_name(1:strfind(sel_name, ')')-1));
            set(ui.chan_name, 'String', p.channel_names{sel_num}, 'Enable', 'on');            
            set(ui.chan_rate, 'Value', min(4,floor((p.channel_rate(sel_num)+4000) / 5000) + 1), 'Enable', 'on');
        end
    end
   
    function calculate_filter_coefficients(index)
        % Calculates the quantized filter coefficients for the given
        % Event index.
        b = [1 0 0]; % Default identity filter.
        a = b;
        fs2 = 10000; % Default sample rate / 2
        source = p.event_source1(index);
        if source == 0
            p.event_source2(index);
        end
        if source > 0
            if p.channel_rate(source) > 0
                fs2 = p.channel_rate(source) / 2;
            end
        end
        if p.event_lower_bandwidth(index) > 0 % Use bandpass filter
            [b, a] = butter(1, [min(fs2-1,p.event_lower_bandwidth(index)) min(fs2-1,p.event_upper_bandwidth(index))] / fs2);
        elseif p.event_upper_bandwidth(index) > p.event_lower_bandwidth(index) % Use lowpass filter
            [b, a] = butter(2, min(fs2-1,p.event_upper_bandwidth(index)) / fs2);
        end
        p.event_filter_a(index,:) = round(a * 2^24); % Quantize filter coefficients
        p.event_filter_b(index,:) = round(b * 2^24);
    end

    function source_index = get_source_index(menu_handle)
        % Return a source channel index. 0 == none, 1 = Ch1, 2 = Ch2, etc
        source_index = 0;
        try
            selection = getnum(menu_handle);  % Selection index
            if selection > 0
                menu_list = get(menu_handle, 'String');
                source_name = menu_list{selection + 1};  % Selected string
                for ilist = 1:length(p.channel_names)
                    if strcmp(source_name, p.channel_names{ilist})
                        source_index = ilist;  % found a match in the channel name list
                        break;
                    end
                end
            end
        catch
        end
    end

    function set_source_index(menu_handle, channel_source_index)
        % Sets the pulldown source menu to a given channel.  If no
        % match is found in the menu, the first menu item is selected.
        menu_index = 0;
        try
            if channel_source_index > 0
                name = p.channel_names{channel_source_index};
                menu_list = get(menu_handle, 'String');
                for ilist = 1:length(menu_list)
                    if strcmp(name, menu_list{ilist})
                        menu_index = ilist; % found channel name in the menu
                        break;
                    end
                end
            end
        catch
        end
        setnum(menu_handle, menu_index - 1);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Callback for Event parameter changed.
    function EventParam(object, event) %#ok<INUSD>
        HighlightSetButton(1);
        index = get(ui.event_list, 'Value');
        list = get(ui.event_list, 'String');
        p.event_method(index) = getnum(ui.event_method);       % 0 = Off, 1 = Window discrim and periodic generator.    
        for i = 1:8
            if p.event_method(i)
                list{i} = ['Event ' num2str(i) ' On'];
            else
                list{i} = ['Event ' num2str(i)];
            end
        end
        set(ui.event_list , 'String', list);
        p.event_source1(index) = get_source_index(ui.event_source1);     % 0 = none, 1 = first analog channel...
        p.event_source2(index) = get_source_index(ui.event_source2); % 0 = none, 1 = first analog channel...
        p.event_transform(index) = getnum(ui.event_transform); % 0 = none, 1 = rectify.
        p.event_lower_bandwidth(index) = getnum(ui.event_lower_bandwidth, [0 ui.fs2-500]);
        p.event_upper_bandwidth(index) = getnum(ui.event_upper_bandwidth, [p.event_lower_bandwidth(index), ui.fs2]);
        calculate_filter_coefficients(index);
        p.event_ratio(index) = getnum(ui.event_ratio, [1 250]);
        p.event_interval_min(index) = round(getnum(ui.event_interval_min, [0 65000]));
        p.event_interval_max(index) = round(getnum(ui.event_interval_max, [0 65000]));
        p.event_interval_exp(index) = 0.1 * floor(10 * getnum(ui.event_interval_exp, [0 500]));
        p.event_threshold(index) = getnum(ui.event_threshold, [-6300 6300]);
        p.event_window1_delay(index) = getnum(ui.event_window1_delay, [0 1000]);
        p.event_window1_max(index) = getnum(ui.event_window1_max, [-6300 6300]);
        p.event_window1_min(index) = getnum(ui.event_window1_min, [-6300 6300]);
        p.event_window2_delay(index) = getnum(ui.event_window2_delay, [0 1000]);
        p.event_window2_max(index) =  getnum(ui.event_window2_max, [-6300 6300]);
        p.event_window2_min(index) = getnum(ui.event_window2_min, [-6300 6300]);
        p.event_refractory_period(index) = round(getnum(ui.event_refractory_period, [0 65000])); % in milliseconds
        p.event_refractory_count(index) = floor(getnum(ui.event_refractory_count, [0 125]));     % If more than 125 are needed, MCU firmware needs to modified.
        p.event_refractory_window(index) = round(getnum(ui.event_refractory_window, [0 65000])); % in milliseconds
        if ismember(p.event_source1(index), ui.channel_Accel_index) || ismember(p.event_source2(index), ui.channel_Accel_index)
            % Accelerometer channels cannot be filtered and only one can be selected.
            p.event_lower_bandwidth(index) = 0;
            p.event_upper_bandwidth(index) = 0;
            if object == ui.event_source1
                p.event_source2(index) = 0;
            elseif object == ui.event_source2
                p.event_source1(index) = 0;
            end
        end
        EventSelect(0,0);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Callback for change to currently selected Event.
    function EventSelect(object, event) %#ok<INUSD>
        if (object == ui.event_list)
            ui.sweep_ok(:) = 0; % Clear sweeps if event_list object was changed.
        end
        index = get(ui.event_list, 'Value');
        list = get(ui.event_list, 'String');
        for i = 1:8
            if p.event_method(i)
                list{i} = ['Event ' num2str(i) ' On'];
            else
                list{i} = ['Event ' num2str(i)];
            end
        end
        set(ui.event_list , 'String', list);
        setnum(ui.event_method, p.event_method(index));
        set_source_index(ui.event_source1, p.event_source1(index));
        set_source_index(ui.event_source2, p.event_source2(index));
        setnum(ui.event_transform, p.event_transform(index));
        setnum(ui.event_lower_bandwidth, p.event_lower_bandwidth(index));
        setnum(ui.event_upper_bandwidth, p.event_upper_bandwidth(index));
        setnum(ui.event_ratio, p.event_ratio(index));
        setnum(ui.event_interval_min, p.event_interval_min(index));
        setnum(ui.event_interval_max, p.event_interval_max(index));
        setnum(ui.event_interval_exp, p.event_interval_exp(index));
        setnum(ui.event_threshold, p.event_threshold(index));
        setnum(ui.event_window1_delay, p.event_window1_delay(index));
        setnum(ui.event_window1_max, p.event_window1_max(index));
        setnum(ui.event_window1_min, p.event_window1_min(index));
        setnum(ui.event_window2_delay, p.event_window2_delay(index));
        setnum(ui.event_window2_max, p.event_window2_max(index));
        setnum(ui.event_window2_min, p.event_window2_min(index));
        setnum(ui.event_refractory_period, p.event_refractory_period(index));
        setnum(ui.event_refractory_count, p.event_refractory_count(index));
        setnum(ui.event_refractory_window, p.event_refractory_window(index));
        if ismember(p.event_source1(index), ui.channel_Accel_index) || ismember(p.event_source2(index), ui.channel_Accel_index)
            if (p.event_source2(index) + p.event_source1(index)) == ui.channel_Accel_index(4)
                set(ui.event_unit_text, 'String', 'C');    % Temperature channel is in Celcius
            else
                set(ui.event_unit_text, 'String', 'mgs');  % Accelerometer scale is 1 gravity / 1000
            end
        else
            if ismember(p.event_source1(index), ui.channel_Aux_index) || ismember(p.event_source2(index), ui.channel_Aux_index)
                set(ui.event_unit_text, 'String', 'mV');     % Millivolts units for Aux channels
            else
                set(ui.event_unit_text, 'String', 'uV');     % Microvolts units for biophysical channels.
            end
        end
        update_plots();
    end

    %%%
    % Callback for change to selected condition.
    function ConditionSelect(object, event) %#ok<INUSD>
        index = getnum(ui.cond_list) + 1; % Index for currently selected condition
        if p.stim_cond_duration(index) <= 0
            set(ui.cond_duration, 'String', 'Unlimited');
            p.stim_cond_duration(index) = 0;
        else
            setnum(ui.cond_duration, p.stim_cond_duration(index));
        end
        setnum(ui.cond_repeat_cond, p.stim_cond_repeat_condition(index));
        setnum(ui.cond_repeat_count, p.stim_cond_repetitions(index));
        if (p.stim_cond_repeat_condition(index) == index) || (p.stim_cond_repeat_condition(index) >= 4)
            set([ui.cond_repeat_count get(ui.cond_repeat_count, 'UserData')], 'Visible', 'off');
        else
            set([ui.cond_repeat_count get(ui.cond_repeat_count, 'UserData')], 'Visible', 'on');
        end
        setnum(ui.cond_stim_mode, p.stim_cond_stim_mode(index));
        setnum(ui.start_cond, p.stim_cond_start);
        setnum(ui.stim_artifact_time, p.stim_artifact_time);
        setnum(ui.stim_aux3_input, p.stim_aux3_input);
        if p.stim_power > 1
            setnum(ui.stim_power, 2); % Select 60V compliance from pull down menu
        else
            setnum(ui.stim_power, max(0, p.stim_power)); % Off or 10V compliance
        end
        setnum(ui.stim_fast_settle, p.stim_fast_settle);
        for i = 1:p.stim_channel_count
            setnum(ui.stim_level(i), p.stim_level(index, i));
            setnum(ui.stim_width(i), p.stim_width(index, i));
            setnum(ui.stim_delay(i), p.stim_delay(index, i));
            setnum(ui.stim_pulses(i), p.stim_pulses(index, i));
            setnum(ui.stim_isi(i), p.stim_isi(index, i));
            setnum(ui.stim_level2(i), p.stim_level2(index, i));
            setnum(ui.stim_width2(i), p.stim_width2(index, i));
            setnum(ui.stim_cathode(i), p.stim_cathodic_first(index, i));
            setnum(ui.stim_anode(i), p.stim_anodic_first(index, i));
            setnum(ui.stim_limit_count(i), p.stim_limit_count(index, i));
            setnum(ui.stim_limit_window(i), p.stim_limit_window(index, i));
            setnum(ui.stim_limit_timeout(i), p.stim_limit_timeout(index, i));
        end
        ConditionParam(0,0);
    end

    %%%
    % Callback for change to a stimulus condition parameter.
    function ConditionParam(object, event) %#ok<INUSD>
        HighlightSetButton(1);
        index = getnum(ui.cond_list) + 1;
        if strcmp(get(ui.cond_duration, 'String'), 'Unlimited')
            p.stim_cond_duration(index) = 0;
        else
            p.stim_cond_duration(index) = getnum(ui.cond_duration, [0 65535]);
            if p.stim_cond_duration(index) == 0
                set(ui.cond_duration, 'String', 'Unlimited');
            end
        end
        p.stim_cond_repeat_condition(index) = getnum(ui.cond_repeat_cond);
        p.stim_cond_repetitions(index) = getnum(ui.cond_repeat_count, [0 65500]);
        if (p.stim_cond_repeat_condition(index) == index) || (p.stim_cond_repeat_condition(index) >= 8)
            set([ui.cond_repeat_count get(ui.cond_repeat_count, 'UserData')], 'Visible', 'off');
        else
            set([ui.cond_repeat_count get(ui.cond_repeat_count, 'UserData')], 'Visible', 'on');
        end
        p.stim_cond_stim_mode(index) = getnum(ui.cond_stim_mode);
        p.stim_cond_start = getnum(ui.start_cond);
        p.stim_artifact_time = getnum(ui.stim_artifact_time, [0 400]);
        p.stim_aux3_input = getnum(ui.stim_aux3_input, [0 8]);
        p.stim_fast_settle = getnum(ui.stim_fast_settle);
        p.stim_power = max(0, getnum(ui.stim_power));
        if (p.stim_power > 1)
            p.stim_power = 6; % Code for 60 V compliance.  Currently there are 'Off', '10 Volt', and '60 V' compliance selections in the pull down menu.
        end
        phase2flag = get(ui.phase2_checkbox, 'Value');
        if phase2flag
            set(ui.phase2_uA_text, 'Visible', 'on');
            set(ui.phase2_width_text, 'Visible', 'on');
        else
            set(ui.phase2_uA_text, 'Visible', 'off');
            set(ui.phase2_width_text, 'Visible', 'off');
        end
        for i = 1:p.stim_channel_count
            p.stim_level(index, i) = getnum(ui.stim_level(i), [-10000 10000]);
            p.stim_width(index, i) = getnum(ui.stim_width(i), [0 10]);
            p.stim_delay(index, i) = getnum(ui.stim_delay(i), [0 6500]);
            p.stim_pulses(index, i) = getnum(ui.stim_pulses(i), [0 250]);
            p.stim_isi(index, i) = getnum(ui.stim_isi(i), [1 6550]);
            p.stim_cathodic_first(index, i) = getnum(ui.stim_cathode(i));
            p.stim_anodic_first(index, i) = getnum(ui.stim_anode(i));
            if (p.stim_cathodic_first(index, i) > 0) && (p.stim_cathodic_first(index, i) ==  p.stim_anodic_first(index, i))
                replacement = [2 1 4 3 6 5];  % Don't allow cathod == anode. Instead pick the matching differential return.
                if object == ui.stim_cathode(i)
                    p.stim_anodic_first(index, i) = replacement(p.stim_anodic_first(index, i));
                    setnum(ui.stim_anode(i), p.stim_anodic_first(index, i));
                else
                    p.stim_cathodic_first(index, i) = replacement(p.stim_cathodic_first(index, i));
                    setnum(ui.stim_cathode(i), p.stim_cathodic_first(index, i));
                end
            end
            p.stim_limit_count(index, i) = floor(getnum(ui.stim_limit_count(i), [0 65000]));
            p.stim_limit_window(index, i) = floor(10 * getnum(ui.stim_limit_window(i), [0 6500])) / 10;
            p.stim_limit_timeout(index, i) = floor(getnum(ui.stim_limit_timeout(i), [0 65000]));
            if phase2flag && (i <= p.stim_cond_stim_mode(index))
                set(ui.stim_level2(i), 'Visible', 'on')
                set(ui.stim_width2(i), 'Visible', 'on')
                p.stim_level2(index, i) = getnum(ui.stim_level2(i), [-10000 10000]);
                p.stim_width2(index, i) = getnum(ui.stim_width2(i), [0 10]);
            else
                set(ui.stim_level2(i), 'Visible', 'off')
                set(ui.stim_width2(i), 'Visible', 'off')
                p.stim_level2(index, i) = getnum(ui.stim_level(i), [-10000 10000]);
                p.stim_width2(index, i) = getnum(ui.stim_width(i), [0 10]);
            end
            if i <= p.stim_cond_stim_mode(index)
                visible = 'on';
            else
                visible = 'off';
            end
            set([ui.stim_level(i) ui.stim_width(i) ui.stim_delay(i) ...
                ui.stim_pulses(i) ui.stim_isi(i) ui.stim_limit_count(i) ...
                ui.stim_limit_window(i) ui.stim_limit_timeout(i)  ...
                ui.stim_cathode(i) ui.stim_anode(i)], 'Visible', visible);
            set(get(ui.stim_delay(i), 'UserData'), 'Visible', visible);
        end
    end

    %%%
    % Called when an item in the Requirement List is selected
    function ReqSelect(object, event) %#ok<INUSD>
        % Get index of selected item
        index = get(ui.req_list, 'Value');
        if isempty(index)
            index = 1;
        else
            index = index(1);
        end

        % Set up the Condition menu with checked Conditions.
        condstrs = get(ui.req_cond, 'String');
        for i=1:p.stim_cond_count
            str = condstrs{i+1};
            if bitand(p.stim_req_active_cond(index), bitshift(1,(i-1))) == 0
                str(1) = ' ';
            else
                str(1) = 'x';
            end
            condstrs{i+1} = str;
        end
        set(ui.req_cond, 'String', condstrs);

        % Set up the Stim menu with checked Stim Channels.
        stimstrs = get(ui.req_output_chans, 'String');
        for i=1:p.stim_channel_count
            str = stimstrs{i+1};
            if bitand(p.stim_req_active_stim(index), bitshift(1,(i-1))) == 0
                str(1) = ' ';
            else
                str(1) = 'x';
            end
            stimstrs{i+1} = str;
        end
        set(ui.req_output_chans, 'String', stimstrs);
        
        % Copy settings to corresponding UI controls.
        setnum(ui.req_cond, 0);
        setnum(ui.req_output_chans, 0);
        setnum(ui.req_source, p.stim_req_source(index));
        setnum(ui.req_type, p.stim_req_type(index));
        setnum(ui.req_n1, p.stim_req_n1(index));
        setnum(ui.req_t1, p.stim_req_t1(index));
        ReqParam(0, 0); 
    end
 
    function ReqParam(object, event) %#ok<INUSD>
        HighlightSetButton(1);
        index = get(ui.req_list, 'Value');
        if isempty(index)
            return;
        end
        index = index(1);

        p.stim_req_source(index) = getnum(ui.req_source);
        p.stim_req_type(index) = getnum(ui.req_type);
        p.stim_req_n1(index) = getnum(ui.req_n1, [0, 125]);  % If more than 125 are needed, the MCU firmware needs to be modified.
        p.stim_req_t1(index) = getnum(ui.req_t1, [0 65000]); % In milliseconds.
               
        reqstr = [num2str(index) ') '];
        sel = getnum(ui.req_cond);
        condstrs = get(ui.req_cond, 'String');
        stimstrs = get(ui.req_output_chans, 'String');
        activestr = '';
        flags = 0;
        for i=1:p.stim_cond_count
            str = condstrs{i+1};
            if i == sel
                if str(1) == 'x'
                    str(1) = ' ';
                else
                    str(1) = 'x';
                end
                condstrs{i+1} = str;
            end
            if str(1) == 'x'
                activestr = [activestr num2str(i) ',']; %#ok<AGROW>
                flags = bitor(flags, bitshift(1,(i-1)));
            end
        end
        p.stim_req_active_cond(index) = flags;
        if sel > 0
            setnum(ui.req_cond, 0);
        end
        if isempty(activestr)
            activestr = 'None';
            reqstr = [reqstr 'None'];
        else
            activestr = ['Condition ' activestr(1:end-1)];
        end
        condstrs{1} = activestr;
        
        stimstrs{1} = 'None';
        stimstr = '';
        flags = 0;
        sel = getnum(ui.req_output_chans);
        for i=1:p.stim_channel_count
            str = stimstrs{i+1};
            if i == sel
                if str(1) == 'x'
                    str(1) = ' ';
                else
                    str(1) = 'x';
                end
                stimstrs{i+1} = str;
            end
            if str(1) == 'x'
                stimstr = [stimstr num2str(i) ',']; %#ok<AGROW>
                flags = bitor(flags, bitshift(1,(i-1)));
            end
        end
        p.stim_req_active_stim(index) = flags;
        if sel > 0
            setnum(ui.req_output_chans, 0);
        end
        if strcmp(activestr, 'None')
            set([ui.req_source ui.req_type ui.req_n1 ui.req_t1], 'Enable', 'off');
        else
            if isempty(stimstr)
                stimstr = 'None';
            else
                stimstr = ['Stim ' stimstr(1:end-1)];
                stimstrs{1} = stimstr;
            end
            set([ui.req_source ui.req_type], 'Enable', 'on');
            eventstr = get(ui.req_source, 'String');
            eventstr = eventstr{p.stim_req_source(index) + 1};
            typestr = get(ui.req_type, 'String');
            typestr = typestr{p.stim_req_type(index) + 1};
            if ~isempty(strfind(typestr, 'occurs'))
                acceptstr = 'Trigger';
                reqstr = [reqstr activestr ': ' acceptstr ' ' stimstr ' when ' eventstr ' ' typestr];
                set(ui.req_t1, 'Enable', 'off', 'String' , '0');
                set(ui.req_n1, 'Enable', 'off', 'String', '1');
                p.stim_req_n1(index) = 1;
                p.stim_req_t1(index) = 0;
            else
                acceptstr = 'Reject';
                reqstr = [reqstr activestr ': ' acceptstr ' ' stimstr ' when ' eventstr ' ' typestr];
                if typestr(end) == 'N'
                    set([ui.req_n1 ui.req_t1], 'Enable', 'on');
                    reqstr = [reqstr(1:end-1) ' ' num2str(p.stim_req_n1(index)) ' for last ' num2str(p.stim_req_t1(index)) ' ms'];
                else
                    set(ui.req_t1, 'Enable', 'on');
                    set(ui.req_n1, 'Enable', 'off');
                    reqstr = [reqstr ' for last ' num2str(p.stim_req_t1(index)) ' ms'];
                end
            end
            set(get(ui.req_output_chans, 'UserData'), 'String', acceptstr);
        end

        set(ui.req_cond, 'String', condstrs);
        set(ui.req_output_chans, 'String', stimstrs);
        reqstrs = get(ui.req_list, 'String');
        reqstrs{index} = reqstr;
        set(ui.req_list, 'String', reqstrs);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Callback for changes to impedance testing parameters
    function ZCheckParam(object, event) %#ok<INUSL>
        HighlightSetButton(1);
        TabSelect(ui.tab_buttons(4), event); % Force update of tab colors.
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Callback for controls that just need to update the plot
    function UpdateCallback(object, event) %#ok<INUSD>
        if (object == ui.ms_scale)
            ui.sweep_ok(:) = 0;  % Clear sweeps if scale box object was changed.
        end
        update_plots();
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Default Callback to handle tests or objects with no callback function
    function DefaultCallback(object, event) %#ok<INUSD>
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Functions for stuffing bytes into the ui.r register array.
    
    % Set unsigned 8-bit integer
    function UByteReg(offset, value)
        if value >= 0
            val = uint32(value);
        else
            val = uint32(int32(value) + 256);
        end
        ui.r(offset+1) = bitand(val, 255);
    end

    % Set signed 8-bit integer
    function ByteReg(offset, value) %#ok<DEFNU>
        UByteReg(offset, value);
    end

    % Unsigned 8-bit integer
    function value = GetUByteReg(offset) %#ok<DEFNU>
        value = ui.r(offset+1);
    end

    % Signed 8-bit integer
    function value = GetByteReg(offset) %#ok<DEFNU>
        val = ui.r(offset+1);
        if val >= 128
            value = int8(int32(val) - 256);
        else
            value = int8(val);
        end
    end
 
    % Set unsigned 16-bit integer
    function UWordReg(offset, value)
        if value >= 0
            val = uint32(value);
        else
            val = uint32(int32(value) + 65536);
        end
        ui.r(offset+1) = bitand(val, 255);
        ui.r(offset+2) = bitand(bitshift(val, -8), 255);
    end

    % Set signed 16-bit integer clipped to int16 range.
    function WordReg(offset, value)
        if value > 32767
            UWordReg(offset, 32767);
        elseif value < -32768
            UWordReg(offset, -32768);
        else
            UWordReg(offset, value);
        end
    end

    % Unsigned 16-bit integer
    function value = GetUWordReg(offset)
        val = uint32(ui.r(offset+1:offset+2));
        value = bitor(bitshift(val(2), 8), val(1));
    end

    % Signed 16-bit integer
    function value = GetWordReg(offset)
        val = GetUWordReg(offset);
        if val >= 32768
            value = int16(int32(val) - 65536);
        else
            value = int16(val);
        end
    end

    % Set unsigned 32-bit integer
    function ULongReg(offset, value)
        if value >= 0
            val = uint32(value);
        else
            val = uint32(double(value) + 4294967296);
        end
        ui.r(offset+1) = bitand(val, 255);
        ui.r(offset+2) = bitand(bitshift(val, -8), 255);
        ui.r(offset+3) = bitand(bitshift(val, -16), 255);
        ui.r(offset+4) = bitand(bitshift(val, -24), 255);
    end

    % Set signed 32-bit integer
    function LongReg(offset, value)
        ULongReg(offset, value);
    end

    % Unsigned 32-bit integer
    function value = GetULongReg(offset)
        val = uint32(ui.r(offset+1:offset+4));
        value = bitor(bitshift(val(4),24), bitor(bitshift(val(3),16), bitor(bitshift(val(2), 8), val(1))));
    end

    % Signed 32-bit integer
    function value = GetLongReg(offset) %#ok<DEFNU>
        val = GetULongReg(offset);
        if val >= 2147483648
            value = int32(double(val) - 4294967296);
        else
            value = int32(val);
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Convert parameters into a register block to send over the IrDA port.
    function SetParameterBlock()
        ui.r = zeros(ui.registerBlockSize,1, 'uint8');
        UByteReg(0, 's'); % Send parameter block command
        UByteReg(1, ui.packet_marker);
        UByteReg(2, ui.packet_marker);
        UByteReg(3, ui.packet_marker);
        UByteReg(4, p.analog_lower_bandwidth_code);
        UByteReg(5, p.analog_upper_bandwidth_code);
        UByteReg(6, p.dsp_high_pass_code);
        if p.intan_channels == 32
            UByteReg(7, 1); % Intan chip ID for 32 channel version
        else
            UByteReg(7, 0); % Intan chip ID for 16 channel version
        end
        UByteReg(8, getnum(ui.event_list)); % Currently selected event channel.
        UByteReg(9, p.stim_cond_start);     % Starting condition.
        UByteReg(10, sweep_decimation());   % Sweep decimation factor
        UByteReg(11, 0); % Filled in by CPU with sampling resolution (1 = 20k, 4 = 5k)
        UByteReg(12, p.stim_power); % Stimulator voltage compliance level 0=off, 1=10V, ..., 6 = 60 Volts
        if ismember(p.stim_fast_settle, 0:10)
            % Fast settle values 0=off, 1:10 are for on each pulse settings
            UByteReg(13, p.stim_fast_settle); % Intan Fast Settle duration 0=off, 1..16 = 200..3200ms on each pulse.
        else
            % Optional fast settle values 11:20 are for post pulse settings
            UByteReg(13, p.stim_fast_settle+6); % Intan Fast Settle duration 17..32 == 200..3200ms pust pulse.
        end
        UByteReg(14, p.stim_aux3_input);       % Intan Aux3 input selection
        UByteReg(15, getnum(ui.z_check_menu)); % 0=Impedance test off, 1=1

        % Session ID and Artifact suppression clock ticks

        offset = 16;
        UWordReg(offset, ui.session_id);
        offset = offset + 2;
        
        UWordReg(offset, floor(40000 * p.stim_artifact_time(1) / 256)); % Artifact suppression convert ms into 40MHz * 256 clock ticks
        offset = offset + 2;

        % Channel rates

        codes = uint32(min(floor(p.channel_rate / 5000), 3)); % Convert sample rate to 0..3 coding

        % Analog Channels 1..16
        ratecode = uint32(0);
        for ichan = 16:-1:1
            ratecode = bitor(bitshift(ratecode, 2), codes(ichan));
        end
        ULongReg(offset, ratecode); % Rate encoding for first 16 A/D channels.
        offset = offset + 4;

        % Analog Channels 17..32 if using 32 channel Intan chip
        ratecode = uint32(0); % Default to recording off.
        if p.intan_channels == 32
            for ichan = 32:-1:17
                ratecode = bitor(bitshift(ratecode, 2), codes(ichan));
            end
        end
        ULongReg(offset, ratecode); % Rate encoding for second set of 16 A/D channels.
        offset = offset + 4;
        
        % Aux Channels.
        ratecode = uint32(0);
        for ichan = 35:-1:33
            ratecode = bitor(bitshift(ratecode, 2), codes(ichan));
        end
        ULongReg(offset, ratecode); % Rate encoding for AUX channels.
        offset = offset + 4;

        % Discriminator paramters

        start_offset = offset;
        for ichan = 1:p.event_count % Event channel index.
            chan = ichan - 1; % Zero based offset.
            offset = start_offset;
            
            calculate_filter_coefficients(ichan);
            
            rate = 20000;  % Sample rate of filter is based on source1 first, source2 second.
            if (p.event_source1(ichan) > 0)
                rate = p.channel_rate(p.event_source1(ichan));
            elseif (p.event_source2(ichan) > 0)
                rate = p.channel_rate(p.event_source2(ichan));
            end
            
            UByteReg(offset + chan, p.event_method(ichan));
            offset = offset + p.event_count;
            
            UByteReg(offset + chan, p.event_transform(ichan));
            offset = offset + p.event_count;
            
            UByteReg(offset + chan, p.event_source1(ichan));
            offset = offset + p.event_count;
            
            UByteReg(offset + chan, p.event_source2(ichan));
            offset = offset + p.event_count;
            
            UByteReg(offset + chan, p.event_ratio(ichan));
            offset = offset + p.event_count;
            
            UWordReg(offset + 2*chan, p.event_lower_bandwidth(ichan));
            offset = offset + 2 * p.event_count;
            
            UWordReg(offset + 2*chan, p.event_upper_bandwidth(ichan));
            offset = offset + 2 * p.event_count;
            
            UWordReg(offset + 2*chan, p.event_interval_min(ichan));        % In milliseconds
            offset = offset + 2 * p.event_count;
            
            UWordReg(offset + 2*chan, p.event_interval_max(ichan));
            offset = offset + 2 * p.event_count;
            
            UWordReg(offset + 2*chan, 10 * p.event_interval_exp(ichan));  % 10x to give extra decimal place
            offset = offset + 2 * p.event_count;
            
            convert2adu = ui.uv2adu;  % Normal conversion to Analog Digitizing Units for Intan chip.
            convertExp = 1;           % Exponent on conversion needed for accelerometer magnitude channel
            source = max(p.event_source1(ichan), p.event_source2(ichan));
            if ismember(source, ui.channel_Accel_index)
                convert2adu = ui.accel2adu;  % Accelerometer channels have different scale
                if source == ui.channel_Accel_index(5)
                    convertExp = 2;   % Accelerometer magnitude must be squared and divided so that (1000 * 1000) = 2^13.
                    convert2adu = (2^13) / (1000^2);  % mag = (x*x + y*y + z*z) / 8, where x,y,z are 10-bit representation of the 16-bit X,Y,Z values.
                elseif source == ui.channel_Accel_index(4)
                    convert2adu = ui.accelTemp2adu;
                end
            end
            
            WordReg(offset + 2*chan, convert2adu * sign(p.event_threshold(ichan)) * abs(p.event_threshold(ichan))^convertExp);
            offset = offset + 2 * p.event_count;
            
            UWordReg(offset + 2*chan, (rate / 1000) * p.event_window1_delay(ichan)); % Convert ms to samples at the filter sample rate.
            offset = offset + 2 * p.event_count;
            
            WordReg(offset + 2*chan, convert2adu * sign(p.event_window1_min(ichan)) * abs(p.event_window1_min(ichan))^convertExp);
            offset = offset + 2 * p.event_count;
            
            WordReg(offset + 2*chan, convert2adu * sign(p.event_window1_max(ichan)) * abs(p.event_window1_max(ichan))^convertExp);
            offset = offset + 2 * p.event_count;
            
            UWordReg(offset + 2*chan, (rate / 1000) * p.event_window2_delay(ichan)); % Convert ms to samples at the filter sample rate.
            offset = offset + 2 * p.event_count;
            
            WordReg(offset + 2*chan, convert2adu * sign(p.event_window2_min(ichan)) * abs(p.event_window2_min(ichan))^convertExp);
            offset = offset + 2 * p.event_count;

            WordReg(offset + 2*chan, convert2adu * sign(p.event_window2_max(ichan)) * abs(p.event_window2_max(ichan))^convertExp);
            offset = offset + 2 * p.event_count;
            
            WordReg(offset + 2*chan, floor(p.event_refractory_period(ichan))); % Refractory period in ms
            offset = offset + 2 * p.event_count;
            
            WordReg(offset + 2*chan, floor(p.event_refractory_count(ichan))); % Number of events in window
            offset = offset + 2 * p.event_count;
            
            WordReg(offset + 2*chan, floor(p.event_refractory_window(ichan))); % Window size in ms
            offset = offset + 2 * p.event_count;
            
            LongReg(offset + 4*chan, p.event_filter_a(ichan, 2));
            offset = offset + 4 * p.event_count;
            
            LongReg(offset + 4*chan, p.event_filter_a(ichan, 3));
            offset = offset + 4 * p.event_count;
            
            gain = 1; % Non-unity filter gain can be used to compensate for overflow or FPGA oddities.
            if (p.event_source1(ichan) > 0) && (p.event_source2(ichan) > 0)
                % If both channel sources are used, then the FPGA divides
                % by two so that (ChAVal - ChBVal)/2 always fits into 16-bits.
                gain = 2 * gain;  % Compensate by boosting the filter gain.
            end
            if ismember(p.event_source1(ichan), ui.channel_Aux_index) || ismember(p.event_source2(ichan), ui.channel_Aux_index)
                % Aux channels normally return 2.4 Volts = 32768 ADU
                % through the discriminator filters.  But we will scale
                % with the filter gain so 2.4 Volts is 2400 in the GUI.
                gain = gain * 2400 / (ui.adu2uv * 32768);
            end
            
            LongReg(offset + 4*chan, gain * p.event_filter_b(ichan, 1));
            offset = offset + 4 * p.event_count;
            
            LongReg(offset + 4*chan, gain * p.event_filter_b(ichan, 2));
            offset = offset + 4 * p.event_count;
            
            LongReg(offset + 4*chan, gain * p.event_filter_b(ichan, 3));
            offset = offset + 4 * p.event_count;
        end

        % Conditon and Stimulus parameters
        
        start_offset = offset;
        for icond = 1:p.stim_cond_count % Condition index.
            cond = icond - 1; % Zero based offset.
            offset = start_offset;
            
            UWordReg(offset + 2*cond, p.stim_cond_duration(icond));
            offset = offset + 2 * p.stim_cond_count;
            
            UByteReg(offset + cond, p.stim_cond_stim_mode(icond));
            offset = offset + p.stim_cond_count;

            UByteReg(offset + cond, p.stim_cond_repeat_condition(icond));
            offset = offset + p.stim_cond_count;
            
            UWordReg(offset + 2*cond, p.stim_cond_repetitions(icond));
            offset = offset + 2 * p.stim_cond_count;

            for istim = 1:p.stim_channel_count
                UWordReg(offset + 2*cond, 10 * p.stim_delay(icond, istim)); % Stored in tenths of a millisecond.
                offset = offset + 2 * p.stim_cond_count;
            end

            for istim = 1:p.stim_channel_count
                WordReg(offset + 2*cond, p.stim_level(icond, istim));
                offset = offset + 2 * p.stim_cond_count;
            end

            for istim = 1:p.stim_channel_count
                UWordReg(offset + 2*cond, 1000 * p.stim_width(icond, istim)); % Stored in microseconds.
                offset = offset + 2 * p.stim_cond_count;
            end
            
            for istim = 1:p.stim_channel_count
                WordReg(offset + 2*cond, p.stim_level2(icond, istim)); % 2nd phase stimulation parameters
                offset = offset + 2 * p.stim_cond_count;
            end

            for istim = 1:p.stim_channel_count
                UWordReg(offset + 2*cond, 1000 * p.stim_width2(icond, istim)); % Stored in microseconds.
                offset = offset + 2 * p.stim_cond_count;
            end

            for istim = 1:p.stim_channel_count
                UWordReg(offset + 2*cond, p.stim_pulses(icond, istim));
                offset = offset + 2 * p.stim_cond_count;
            end

            for istim = 1:p.stim_channel_count
                UWordReg(offset + 2*cond, 10 * p.stim_isi(icond, istim)); % Stored in tenths of a millisecond.
                offset = offset + 2 * p.stim_cond_count;
            end

            % Stim cathodic/anodic first settings

            for istim = 1:p.stim_channel_count
                UByteReg(offset + cond, p.stim_cathodic_first(icond, istim));
                offset = offset + p.stim_cond_count;
            end
            for istim = 1:p.stim_channel_count
                UByteReg(offset + cond, p.stim_anodic_first(icond, istim));
                offset = offset + p.stim_cond_count;
            end

            % Stimulation limit parameters

            for istim = 1:p.stim_channel_count
                UWordReg(offset + 2*cond, floor(p.stim_limit_timeout(icond, istim)));     % In seconds
                offset = offset + 2 * p.stim_cond_count;
            end

            for istim = 1:p.stim_channel_count
                UWordReg(offset + 2*cond, floor(p.stim_limit_count(icond, istim)));
                offset = offset + 2 * p.stim_cond_count;
            end

            for istim = 1:p.stim_channel_count
                UWordReg(offset + 2*cond, floor(p.stim_limit_window(icond, istim) * 10)); % In 100 ms increments
                offset = offset + 2 * p.stim_cond_count;
            end
        end

        % Stimulation Requirement parameters
        
        start_offset = offset;
        for ireq = 1:16 % Requirement index
            req = ireq - 1;  % Zero based offset.
            offset = start_offset;
            
            UByteReg(offset + req, p.stim_req_active_cond(ireq));
            offset = offset + p.stim_req_count;
            
            UByteReg(offset + req, p.stim_req_active_stim(ireq));
            offset = offset + p.stim_req_count;
            
            UByteReg(offset + req, p.stim_req_source(ireq));
            offset = offset + p.stim_req_count;
            
            UByteReg(offset + req, p.stim_req_type(ireq));
            offset = offset + p.stim_req_count;
            
            UWordReg(offset + 2*req, p.stim_req_n1(ireq));
            offset = offset + 2 * p.stim_req_count;
            
            UWordReg(offset + 2*req, p.stim_req_t1(ireq));
            offset = offset + 2 * p.stim_req_count; %#ok<NASGU>
        end

        % Checksum for block (excluding header and checksum)
        ULongReg(ui.registerBlockSize-4, sum(ui.r(5:ui.registerBlockSize-4)));
    end


end % end main function (because nested function was used, 
    % an end is required on all functions/