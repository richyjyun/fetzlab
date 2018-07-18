% Window Discrimination visualization for Synapse
% Irene Rembado, Larry Shupe - University of Washington 2017

close all; clear all; clc;

addpath('C:\TDT\OpenEx\Examples\TTankX_Example\Matlab')
addpath('C:\TDT\Synapse\SynapseAPI\Matlab')

%% Initialize parameters. 
% Single sweep reading from TDT centered on Threshold crossing 
%%%%  They have to be the same as in StStore1 object in Synapse %%%%
Duration = 200; % msec
PreCapture = 10; % msec
fs = 24414; %sampling rate

% based on Strobe Store gizmo
EVENT = 'StS1';
StimID = 'Thd1';
Thr = 'Thdb';
Delay1 = 'Thdc';
W1max = 'Thdd';
W1min = 'Thde';
Delay2 = 'Thdf';
W2max = 'Thdg';
W2min = 'Thdh';

% Init parameters for plot
fig_color = [1 1 1];
y_low = -400;% unity (millivolt)
y_high = 400;% unity (millivolt)
Buffer_length = Duration/1000*fs;
Buffer_time = (Duration-PreCapture);
time_zero = round(PreCapture/1000*fs);
zero_array = zeros(1,Buffer_length);

t = SynapseLive();
t.TYPE = {'snips', 'epocs', 'scalars'};
t.VERBOSE = false;

first_pass = true;

%% init Plot
screen_size = [0 0 1500 1024];
fig_title = strcat('Wind Discrim');
figure_handle = figure('Name',fig_title,'MenuBar','none','NumberTitle','off','Color',fig_color,'Position',[screen_size(4)/2 screen_size(4)/4 screen_size(3)/1.5 screen_size(4)/1.5]);
axes_handles = axes('Position',[0.1,0.1,0.7,0.8]);   
plot_handles = plot(zero_array);
hold on
line([time_zero time_zero],[y_low y_high],'color','k')
% sets fontsize
set(axes_handles,'FontSize',12,'XTick',[0, Buffer_length ],'XTickLabel',{['-',num2str(PreCapture)], num2str(Buffer_time) });
% set y axis limit
ylim([y_low,y_high]);
ylabel('mV');
% set default x axis limit
xlim([0,Buffer_length]);
xlabel('msec')

% xlim ms editable text
% data are in b_xlim.String
x_lim_ms_init = Buffer_time;
b_name = num2str(x_lim_ms_init);
b_xlim_ms = uicontrol('Style','edit','Units','normalized','String', b_name ,'Value',0,'FontSize',13,'Position',[0.85,0.75,0.08,0.05]);
text_xlim = uicontrol('Style','text','Units','normalized','String', 'Xscale(ms)','FontSize',13,'Position',[0.85,0.81,0.09,0.03],'BackgroundColor','w');
x_lim_ms_update = false;

% ylim editable text
% data are in b_ylim.String
y_lim_init = y_high;
b_name = num2str(y_lim_init);
b_ylim = uicontrol('Style','edit','Units','normalized','String', b_name ,'Value',0,'FontSize',13,'Position',[0.85,0.6,0.08,0.05]);
text_ylim = uicontrol('Style','text','Units','normalized','String', 'Yscale(mV)','FontSize',13,'Position',[0.85,0.66,0.09,0.03],'BackgroundColor','w');
y_lim_update = false;

% clear button
b_clear = 'Clear';
clear_value_init = 0;
b_clear = uicontrol('Style','ToggleButton','Units','normalized','Value',clear_value_init,'String',b_clear,'FontSize',10,'FontWeight','bold','Position',[0.85,0.1,0.08,0.05]);
clear_update = false;  

while 1
    
    % slow it down a little
    pause(1) % this value can't be too small otherwise it misses data points
    
    % get the most recent data
    t.update;   
   
    % get snippet events and epochs
    r = t.get_data(EVENT);
    stimID = t.get_data(StimID);
    thr = t.get_data(Thr);
    delay1 = t.get_data(Delay1);
    w1max = t.get_data(W1max);
    w1min = t.get_data(W1min);
    delay2 = t.get_data(Delay2);
    w2max = t.get_data(W2max);
    w2min = t.get_data(W2min);
       
    % check xlim ms
    x_lim_ms_now = str2double(get(b_xlim_ms,'String'));
    if x_lim_ms_init ~= x_lim_ms_now
        if (x_lim_ms_now<=Buffer_time) && (x_lim_ms_now>0)
            x_lim_ms_init = x_lim_ms_now;
            x_lim_ms_update = true;
        else
            % button reset
            set(b_xlim_ms,'String', num2str(Buffer_time))
        end
    end
    
    % check ylim
    y_lim_now = str2double(get(b_ylim,'String'));
    if y_lim_init ~= y_lim_now
        y_lim_init = y_lim_now;
        y_lim_update = true;
    end
    % check clear button 
     clear_now = get(b_clear,'Value');
    if clear_now ~= clear_value_init
        clear_value_init = clear_now;
        clear_update = true;
    end
    
    % update plot
    if isstruct(r) && isstruct(thr) 
       if ~isnan(r.data)
            
           % get data from channel 1 (It requires the selection of one channel in Synapse)
           ind = find(r.chan == 1);
           nstrobes = size(r.data(ind,:),1);
           for nsize = 1:nstrobes
               ddd = r.data(nsize,:);

                if first_pass
                    nsweeps = 0;
                    first_pass = false;
                    plot_data = ddd;
                else
                    plot_data = ddd;
                end

                fs = r.fs; % sampling rate
                dey1 = (delay1.data+PreCapture)./1000*fs;
                dey2 = (PreCapture+delay2.data)/1000*fs;

                % plot single sweep data and 
                tdiff = abs(stimID.onset - r.ts(nsize));
                [mindiff, ind] = min(tdiff);
                if stimID.data(ind) == 1 % timestamp of sweep nearest timestamp threshold crossing
                      plot(plot_data,'color','g')
                else
                      plot(plot_data,'color',[0.8 0.8 0.8])
                end
            end
           
            nsweeps = nsweeps + nstrobes;
                      
            thr_line = line([0 Buffer_length], [thr.data(end) thr.data(end)],'color','r');
            line([dey1(end) dey1(end)], [w1min.data(end) w1max.data(end)],'color','r')
            line([dey2(end) dey2(end)], [w2min.data(end) w2max.data(end)],'color','r')
            ttt = sprintf('nsweeps = %d', nsweeps);
            title(ttt)
            
            % update x axis
           if x_lim_ms_update 
                x_lim_ms_update = false;
                x_limit_now = round((x_lim_ms_now+PreCapture)/1000*fs);
                % buffer length workaround
                if ( x_limit_now > Buffer_length )
                    x_limit_now = Buffer_length;
                    x_lim_ms_now = Buffer_time;
                    display('The xlim cannot exceed the buffer length')
                    set(b_xlim_ms,'String',num2str(Buffer_time))
                end
                set(axes_handles,'xlim',[0, x_limit_now]);
                % update xtick
                set(axes_handles,'FontSize',12,'XTick',[0, x_limit_now ],'XTickLabel',{['-',num2str(PreCapture)], num2str(x_lim_ms_now) });
                drawnow;
           end
        % update y axis
        if y_lim_update 
            y_lim_update = false;
            y_low = -y_lim_now;
            y_high = y_lim_now;
            set(axes_handles,'ylim',[y_low,y_high]);
            drawnow;
        end
        % clear the plot if clear button 
        if clear_update || (nsweeps > 2000)
            nsweeps = 0;
            clear_update = false;
            cla(axes_handles)
            plot_handles=plot(zero_array);
            line([time_zero time_zero],[y_low y_high],'color','k');
            drawnow;
        end
                   
            % plot refresh
            drawnow
        end
    end
end


