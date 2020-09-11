% fucntion used in fancyPlot to get indices of when sleep states begin and
% end
function[state_indxOn, state_indxOff] = getStates(codes,state_times,time_ds,gcamp_norm, stage);

%%get state change times for plots
wake_indx = find(codes == stage);

if wake_indx(1)-1 <  0   %AV changed 10th Jan as wasn't getting first sleep episode
    wake_indx = wake_indx(2:end);
end

wake = nonzeros(wake_indx);
wakePlus = wake+1;
%wake = wake + 2; uncomment this line to get transitions FROM a state

%wake = wake(2:end-1);
%sws = find(codes == 2);
wake_times = state_times(wake);

wake_timesPlus = zeros(length(wakePlus),1);
if wakePlus(end)>length(codes)
    wake_timesPlus(1:end-1)=state_times(wakePlus(1:end-1));
    wake_timesPlus(end) = time_ds(end);
else
    wake_timesPlus=state_times(wakePlus);
end

%wake_timesPlus = state_times(wakePlus);
wake_times2 = floor(wake_times);
wake_times2Plus = floor(wake_timesPlus);
wake_times1 = nonzeros(wake_times2);
wake_times1Plus = nonzeros(wake_times2Plus);

state_indxOn = zeros(length(wake_times1),1);
state_indxOff = zeros(length(wake_times1Plus),1);

for i = 1:length(wake_times1)
    state_indxOn(i) = find(time_ds>wake_times1(i),1);
end

for i = 1:length(wake_times1Plus)
    state_indxOff(i) = find(time_ds>wake_times1Plus(i),1);
    if wake_times1Plus(i)>time_ds(end)
        state_indxOff(i) = find(time_ds(end));
    end
end
