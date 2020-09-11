% Data to analyse
mouseName = 'yew1501_fancy';
Mouse = strcat(mouseName,'.mat');
mouse_coeff = strcat('coeff_',Mouse);

% User-defined settings
normTag = 1; % 0 for normal %, 1 for z-score, 2 for percentage of min and max point)
plotFit = 0;  %if 1,  save fitted curve
uvCorrectionTag = 0;  % If 1, peform background subtraction 
Acq_rate = 1; % user input (Hz): choose time resolution of final figure

% Makes a 'results' folder in current directory in which to export
% data/figures
current_folder = pwd;
if ~exist('averages', 'dir')
    mkdir('averages')
end

export_folder = strcat(current_folder,  '\averages\');

%% Load data from mouse
load(Mouse)

%get Ca2+ and UV signal from loaded data
gcamp_signal = Calcium.values;
%uv_signal = Backgrou.values;

% Get sleep stage (state) marker codes and times
state_times = instant_state.times;
codes = instant_state.codes(:,1);

%Extracting timing data from Calcium struct
N = Calcium.length; % Number of datapoints
dt = Calcium.interval; % (Sampling interval from each data point in Calcium struct) 
tmin = Calcium.start; % First Time value (in seconds)
tmax = Calcium.start + (N*dt); % Time at the end of data (in seconds)
t = linspace(tmin, tmax, N); % Equi-distant time vector (for FFT)

% Determine if there is a light stimulus and only analyse up to this point
% if so

lightPulse = find(codes == 7); % find index of light pulse
if size(lightPulse)>0
  lightTag = 1;
else
    lightTag = 0;
end

if lightTag ==1
    lightTime = state_times(lightPulse);
    lightOnIndx = find(t>lightTime,1)-1;
    gcamp_signal = gcamp_signal(1:lightOnIndx);
    uv_signal = uv_signal(1:lightOnIndx);
    t = t(1:lightOnIndx);
else
end

% Determine if there are NaNs in the data
nanCheck = sum(isnan(Calcium.values));
if nanCheck>0
    nanTag = 1;
else
    nanTag = 0;
end

%down sample Calcium and UV signal to 'Acq_rate' by taking the mean every 'ds_factor' samples (in our case)
ds_factor = (1/dt)/Acq_rate;
gcamp_ds = group_z_project_vector(gcamp_signal, ds_factor);    
%uv_ds = group_z_project_vector(uv_signal, ds_factor);
time_ds = group_z_project_vector(t', ds_factor);

%% Fit the Calcium and UV data to 2nd order exponential curves
if nanTag == 1
    gcampValid = ~isnan(gcamp_ds);
    %uvValid = ~isnan(uv_ds);
    f=fit(time_ds(gcampValid)', gcamp_ds(gcampValid)','exp2');   % exp2: Y = a*exp(b*x)+c*exp(d*x)
    %f2=fit(time_ds(uvValid)', uv_ds(uvValid)','exp2');
else
    f=fit(time_ds', gcamp_ds','exp2');   % exp2: Y = a*exp(b*x)+c*exp(d*x)
    %f2=fit(time_ds', uv_ds','exp2');
end

% Plot Calcium and UV signal, together with fitted curves
figure(1);
hold on;
plot(f,'r',time_ds',gcamp_ds','b');
%plot(f2,'r',time_ds',uv_ds','m');
hold off;



% Calculate coefficients of the fitted curves
coeffs1 = coeffvalues(f);  % Calcium coefficients
%coeffs2 = coeffvalues(f2);  %UV coefficients

% Re-create curves using the same timebase as Calcium and UV signals s cna
% perform dF/F calculation
fitData = coeffs1(1)*exp(coeffs1(2)*time_ds) + coeffs1(3)*exp(coeffs1(4)*time_ds);
%fitData2 = coeffs2(1)*exp(coeffs2(2)*time_ds) + coeffs2(3)*exp(coeffs2(4)*time_ds);

% Plot Calcium and UV signal, together with fitted curves
figure(1);
hold on;
calcium_plot = plot(time_ds',gcamp_ds','b');
%calcium_fit = plot(f,'k');
%uv_plot = plot(time_ds',uv_ds','m');
%uv_fit = plot(f2,'r');

if exist(mouse_coeff,'file')
    load(mouse_coeff);
    %fitData = fittedmodel.a*exp(fittedmodel.b*time_ds) + fittedmodel.c*exp(fittedmodel.d*time_ds);
    fitData = transpose(fittedmodel(time_ds));   % 'fittedmodel' includes algotihm and coefficients from curve-fitting session
    calcium_userfit = plot(time_ds',fitData','k--', 'LineWidth',2);
    legend([calcium_plot calcium_userfit],{'Calcium','Calcium Userfit'});
    %legend([calcium_plot calcium_userfit uv_plot uv_fit],{'Calcium','Calcium Userfit', 'UV', 'UV Autofit'});
else
    fitData = coeffs1(1)*exp(coeffs1(2)*time_ds) + coeffs1(3)*exp(coeffs1(4)*time_ds);
    calcium_fit = plot(time_ds',fitData','k');
    legend([calcium_plot calcium_fit],{'Calcium','Calcium Autofit'});
    %legend([calcium_plot calcium_fit uv_plot uv_fit],{'Calcium','Calcium Autofit', 'UV', 'UV Autofit'})
end

hold off;

if plotFit == 1
    saveas(gcf,strcat(export_folder,mouseName,'_fittedData.png'));
end

% Calculate dF/F
if nanTag == 1
    normData = (gcamp_ds(gcampValid) - fitData(gcampValid))./fitData(gcampValid);
    %normData2 = (uv_ds(uvValid) - fitData2(uvValid))./fitData2(uvValid);
else
    normData = (gcamp_ds - fitData)./fitData;
    %normData2 = (uv_ds - fitData2)./fitData2;
end

if uvCorrectionTag == 0
    signalNorm2 = normData;
else 
    signalNorm2 = normData-normData2;
end

if plotFit == 1
    saveas(gcf,strcat(export_folder,mouseName,'_fittedData.png'));
end

%% Normalize data based on normTag
minSignalNorm = min(signalNorm2);
maxSignalNorm = max(signalNorm2);
diffSignalNorm = maxSignalNorm - minSignalNorm;
stdData = nanstd(signalNorm2);
meanData = nanmean(signalNorm2);
zData = (signalNorm2-meanData)/stdData;

if normTag == 0
    signalNorm = signalNorm2*100;
else
    if normTag == 1
        signalNorm = zData;
    else
        if normTag == 2
            signalNorm =(signalNorm2-minSignalNorm)/diffSignalNorm*100;
        end
    end
end

% Plot dF/F for Calcium and UV, together with normalized data
figure(3)
hold on;
if nanTag == 1
    plot(time_ds (gcampValid), normData, 'b')
    %plot(time_ds(uvValid), normData2, 'm')
    signalNormValid = ~isnan(signalNorm);
    plot(time_ds(signalNormValid), signalNorm, 'r')     
else
    plot(time_ds, normData, 'b')
    %plot(time_ds, normData2, 'm')
    plot(time_ds, signalNorm, 'r')
end 
hold off;


%% Get state change times for plots
gcamp_norm = signalNorm;

%% Find code of the light pulse and don't go past this
all_codes = find(codes<10);

for i = 1:length(codes)
    if all_codes(i) > lightPulse
        all_codes = codes(1:i-1);
    else 
        all_codes = codes;
    break
    end
end

state_timesAcq = (instant_state.times- tmin)*Acq_rate; % gives times in Acq_rate starting at 0

if lightTag == 1
    timePost = state_timesAcq(lightPulse);
else
    timePost = state_timesAcq(end);
end
timePre = state_timesAcq(1);

[wakeAvgs,wakeAUCs,wake_timesIntegers,figure3, figure4] = getStateAvgs(timePre,timePost,1,state_timesAcq,all_codes,gcamp_norm);
[nremAvgs,nremAUCs,nrem_timesIntegers,figure3, figure4] = getStateAvgs(timePre,timePost,2,state_timesAcq,all_codes,gcamp_norm);
[remAvgs,remAUCs,rem_timesIntegers,figure3, figure4] = getStateAvgs(timePre,timePost,3,state_timesAcq,all_codes,gcamp_norm);

save(strcat(export_folder,mouseName,'_wake_avgs.mat'),'wakeAvgs');
save(strcat(export_folder,mouseName,'_nrem_avgs.mat'),'nremAvgs');
save(strcat(export_folder,mouseName,'_rem_avgs.mat'),'remAvgs');

save(strcat(export_folder,mouseName,'_wake_AUCs.mat'),'wakeAUCs');
save(strcat(export_folder,mouseName,'_nrem_AUCs.mat'),'nremAUCs');
save(strcat(export_folder,mouseName,'_rem_AUCs.mat'),'remAUCs');