% Plot photometry figure, colour-coding EEG, EMG and calcium signal by state
% AV: 8th April 2020
% 1) This script uses an exponential fit of both the calcium and the UV
% channel to calculate f-f0/f0.
% 4) This scipt works best if you downsample to 1 Hz (seems to get rid of a
% lot of the 'fuzz'
% 5) If lightTag =1. The co-eff bit refers to using a 'human-fitted' 2nd order exponential
% curve to get the data (i.e. to exclude the high peaks seen during the
% light pulse. This is saved as 'mouseName_coeff.mat' file which is loaded
% with the data here.
% 6) This script also has nanTag (=1 if I have had to NaN any values in the
% data previously (e.g. because LED got turned off) and normTag- which
% either leaves the data as is, calculates z score across data or
% calculates percentage using maximum and minimum of data).

%% Nan Commands
%Calcium.values(Calcium.values<0.67) = nan; %%use this command when need to nan values from a data set
%Backgrou.values(Backgrou.values<0.06) = nan;

%% IMPORTANT!!! 
% This script requires these functions in order to work:
% - getStates.m
% - plotStates.m
% Also required are co-effiicents for cruve fitting (use curve-fitting
% toolbox, save f variables as 'mouseName_coeff.mat'
% Remember to run "set(gcf,'renderer','painters');" to export as emf

%% set expt variables
mouseName = uigetfile;
mousestring = split(mouseName,'.');
mousestring = mousestring{1};
mouse_coeff = strcat('coeff_',mousestring,'.mat');

% Dialog boxes
% 1. To determine whether to account for NaN values
yesno = {'No','Yes'};
[nanTag,~] = listdlg('PromptString',{'Are there NaN values in your data?',''},...
'SelectionMode','single','ListString',yesno);

% 2. To select normalization method
norm_method = {'Percent','Z-score','Min-Max'};
[normTag,~] = listdlg('PromptString',{'What normalization method do you want to use?',''},...
'SelectionMode','single','ListString', norm_method); %(1 for percent, 2 for z-score, 3 for percentage of min and max point)

% 3. To determine time resolution of downsampled data
ds_ui_prompt = {'Choose time resolution in Hz for final figure:'};
ds_dlgtitle = 'Input';
ds_ui_dims = [1 35];
ds_definput = {'1'};
ds_answer = inputdlg(ds_ui_prompt,ds_dlgtitle,ds_ui_dims,ds_definput);
Acq_rate = str2double(ds_answer{1});

%% Downsample data and fit curve to calculate deltaF/F

%get Ca2+ and UV signal from loaded data
data = load(mouseName);
gcamp_signal = data.Calcium.values;
uv_signal = data.Backgrou.values;

% Caclulate the downsampling factor from sampling interval
dt = data.Calcium.interval; % Time intervals 
ds_factor = (1/dt)/Acq_rate; % calculate the downsampling factor

%gcamp_signal = signalNorm;
gcamp_ds = group_z_project_vector(gcamp_signal, ds_factor);    %down sample to 'Acq_rate' by taking the mean every 10 samples (in our case)
uv_ds = group_z_project_vector(uv_signal, ds_factor);

%Extracting time data to construct timeline
N = data.Calcium.length; % Number of datapoints
tmin = data.Calcium.start; % First Time value (in seconds)
tmax = data.Calcium.start + (N*dt); % Time at the end of data (in seconds)
t = linspace(tmin, tmax, N); % Equi-distant time vector (for FFT)

time_ds = group_z_project_vector(t', ds_factor);

%% 
if nanTag == 2
    gcampValid = ~isnan(gcamp_ds);
    uvValid = ~isnan(uv_ds);
    f=fit(time_ds(gcampValid)', gcamp_ds(gcampValid)','exp2');   % exp2: Y = a*exp(b*x)+c*exp(d*x)
    f2=fit(time_ds(uvValid)', uv_ds(uvValid)','exp2');
else
    f=fit(time_ds', gcamp_ds','exp2');   % exp2: Y = a*exp(b*x)+c*exp(d*x)
    f2=fit(time_ds', uv_ds','exp2');
end

figure(1);
hold on
calcium_plot = plot(time_ds',gcamp_ds','b');
calcium_fit = plot(f,'k');
uv_plot = plot(time_ds',uv_ds','m');
uv_fit = plot(f2,'r');

coeffs1 = coeffvalues(f);
coeffs2 = coeffvalues(f2);

if exist(mouse_coeff,'file')
    load(mouse_coeff);
    %fitData = fittedmodel.a*exp(fittedmodel.b*time_ds) + fittedmodel.c*exp(fittedmodel.d*time_ds);
    fitData = transpose(fittedmodel(time_ds));   % 'fittedmodel' includes algotihm and coefficients from curve-fitting session
    calcium_userfit = plot(time_ds',fitData','k--', 'LineWidth',2);
    legend([calcium_plot calcium_fit calcium_userfit uv_plot uv_fit],{'Calcium','Calcium Autofit','Calcium Userfit', 'UV', 'UV Autofit'}); 
else
    [curveTag,tf] = listdlg('PromptString',{'Are you satisfied with the fit to gcamp_ds?',''},...
    'SelectionMode','single','ListString',yesno);
    if curveTag == 1
        f = msgbox(['TO CHANGE THE FITTED LINE OF FIGURE 1, OPEN THE Curve Fitting App NOW AND DO THE FOLLOWING:' newline ...
        '1. Plot "time_ds" as X data and "gcamp_ds" as Y data.' newline...
        '2. Select "Exponential" with "2" terms as fit parameters'...
        '(2nd order exponentials are generally a good fit'...
        ' for bleaching photometry signals, though you can try other fits)' newline...
        '3. Use the exclude outlier button to interactively exclude outliers'...
        ' from being used to calculate the fitted line.' newline...
        '4. Go to the "Fit" menu on the Menu bar and "Save to Workspace".' newline...
        '5. Save the variable, "fittedmodel" in the workspace as "coeff_filename.mat"']);

        return
    else
        fitData = coeffs1(1)*exp(coeffs1(2)*time_ds) + coeffs1(3)*exp(coeffs1(4)*time_ds);     
        calcium_fit = plot(time_ds',fitData','k');
        legend([calcium_plot calcium_fit uv_plot uv_fit],{'Calcium','Calcium Autofit', 'UV', 'UV Autofit'})  
    end
end

hold off

fitData2 = coeffs2(1)*exp(coeffs2(2)*time_ds) + coeffs2(3)*exp(coeffs2(4)*time_ds);

if nanTag == 2
    normData = (gcamp_ds(gcampValid) - fitData(gcampValid))./fitData(gcampValid);
    normData2 = (uv_ds(uvValid) - fitData2(gcampValid)./fitData2(gcampValid));
else
    normData = (gcamp_ds - fitData)./fitData;
    normData2 = (uv_ds - fitData2)./fitData2;
end

signalNorm2 = normData;
%signalNorm2 = normData-normData2; %uncomment this line to get background channel subtraction
minSignalNorm = min(signalNorm2);
maxSignalNorm = max(signalNorm2);
diffSignalNorm = maxSignalNorm - minSignalNorm;
stdData = nanstd(signalNorm2);
meanData = nanmean(signalNorm2);
zData = (signalNorm2-meanData)/stdData;

if normTag == 1
    signalNorm = signalNorm2*100;
else
    if normTag == 2
        signalNorm = zData;
    else
        if normTag == 3
            signalNorm =(signalNorm2-minSignalNorm)/diffSignalNorm*100;
        end
    end
end

gcamp_norm = signalNorm;

%% Get times of state changes and plot gcamp signal colour coded
%%get state change times for plots
codes = data.instant_state.codes(:,1);
state_times = data.instant_state.times; % gives times in seconds

wake_code = 1;
sws_code = 2;
rem_code = 3;

[wake_indxOn, wake_indxOff] = getStates(codes,state_times,time_ds,gcamp_norm,wake_code);
[sws_indxOn,sws_indxOff] = getStates(codes,state_times,time_ds,gcamp_norm,sws_code);
[rem_indxOn, rem_indxOff] = getStates(codes,state_times,time_ds,gcamp_norm,rem_code);

% Make a time x-axis
if t(end)< 3600
    minuteIncrement = 5;
else
    minuteIncrement = 30;
end

timeMin = state_times(1)-1;
timeMax = t(end);
timeInterval = 60*minuteIncrement;
timeMinutesRange = timeMin: timeInterval: timeMax;
timeMinutesLabels = string((timeMinutesRange -timeMin)/60);

fig = figure(2);
ax1 = subplot(4,1,1);

All_wake_trials = cell(length(wake_indxOn));
All_wake_times = cell(length(wake_indxOn));
for i = 1:length(wake_indxOn)
    hold on
    All_wake_trials{i} = gcamp_norm(wake_indxOn(i):wake_indxOff(i)); % deltaF /F(t)  = F(t) - mean(F)/mean(F)
    All_wake_times{i} = time_ds(wake_indxOn(i):wake_indxOff(i));
    wake_plot = plot(All_wake_times{i},All_wake_trials{i},'b','LineWidth',0.5);
end
 
All_sws_trials = cell(length(sws_indxOn));
All_sws_times = cell(length(sws_indxOn));
for i = 1:length(sws_indxOn)
    All_sws_trials{i} = gcamp_norm(sws_indxOn(i):sws_indxOff(i)); % deltaF /F(t)  = F(t) - mean(F)/mean(F)
    All_sws_times{i} = time_ds(sws_indxOn(i):sws_indxOff(i));
    sws_plot = plot(All_sws_times{i},All_sws_trials{i},'g','LineWidth',0.5);
end

All_rem_trials = cell(length(rem_indxOn));
All_rem_times = cell(length(rem_indxOn));

for i = 1:length(rem_indxOn)
    All_rem_trials{i} = gcamp_norm(rem_indxOn(i):rem_indxOff(i)); % deltaF /F(t)  = F(t) - mean(F)/mean(F)
    All_rem_times{i} = time_ds(rem_indxOn(i):rem_indxOff(i));
    rem_plot = plot(All_rem_times{i},All_rem_trials{i},'c','LineWidth',0.5);
end

ax1.XLim = ([state_times(1),max(time_ds)]);
ax1.YLim = ([-2 6]);
ax1.XTick = timeMinutesRange;
ax1.XTickLabel = [];
ax1.YLabel.String = '\deltaF/F (z-score)';
hold off

eeg_start = find(t > state_times(1),1)-1;
eeg_trunc = data.EEG1.values(eeg_start:end);

ax2 = subplot(4,1,2);
spectrogram(eeg_trunc,1024,120,1024,500,'yaxis');
cbar = colorbar;
cbar.Visible = 'on';
ax2.Position=[0.1300 0.5482 0.7750 0.1577];
ax2.YLabel.String = 'FFT (Hz)';
ax2.YLim = ([0,20]);
ax2.CLim = [-50 1];
ax2.XTick = timeMinutesRange;
ax2.XTickLabel = [];
ax2.XLabel = [];

eeg = data.EEG1.values;
[wake_indxOnAll, wake_indxOffAll] = getStates(codes,state_times,t,eeg,wake_code);
[sws_indxOnAll,sws_indxOffAll] = getStates(codes,state_times,t,eeg,sws_code);
[rem_indxOnAll, rem_indxOffAll] = getStates(codes,state_times,t,eeg,rem_code);
ax3 = subplot(4,1,3);

eeg_start = find(t > state_times(1),1)-1;
emg_time = t(eeg_start:end);

wake_eeg = plotStates(wake_indxOnAll, wake_indxOffAll,t,eeg,ax3, 'b');
sws_eeg = plotStates(sws_indxOnAll, sws_indxOffAll,t,eeg,ax3, 'g');
rem_eeg = plotStates(rem_indxOnAll,rem_indxOffAll,t,eeg,ax3, 'c');

ax3.XLim = ([state_times(1),max(time_ds)]);
ax3.XTick = timeMinutesRange;
ax3.XTickLabel = [];
ax3.YLim = ([-5,5]);
ax3.YLabel.String = 'EEG (V)';

ax4 = subplot(4,1,4);
emg = data.EMG1.values;
wake_emg = plotStates(wake_indxOnAll, wake_indxOffAll,t,emg,ax4, 'b');
sws_emg = plotStates(sws_indxOnAll, sws_indxOffAll,t,emg,ax4, 'g');
rem_emg = plotStates(rem_indxOnAll,rem_indxOffAll,t,emg,ax4, 'c');

ax4.XLim = ([state_times(1),max(time_ds)]);
ax4.YLim = ([-10,10]);
ax4.YLabel.String = 'EMG (V)';
ax4.XTick = timeMinutesRange;
ax4.XTickLabel = timeMinutesLabels;

%set(gcf,'renderer','painters'); To make sure this saves and exports as an
%emf (vector graphics object. WARNING: the file is huge and takes ages to
%import...)
%fig.PaperPosition = [0 0 4.53 2.4];
%print(fig,'-dtiff','-r600');