% function to plot a particular state in a particular color
function[state_trials] = plotStates(state_indxOn, state_indxOff,t,data,ax, color);
for i = 1:length(state_indxOn)
    hold on
    %f1(i) = nanmedian(gcamp_norm(wake_times3(i)-pre_state_time:wake_times3(i)+ post_state_time));
    state_trials{i} = data(state_indxOn(i):state_indxOff(i)); % deltaF /F(t)  = F(t) - mean(F)/mean(F)
    All_wake_times{i} = t(state_indxOn(i):state_indxOff(i));
    state_plot = plot(ax,All_wake_times{i},state_trials{i},color,'LineWidth',0.1);
    %All_f_trials(i,:) = gcamp_norm((wake_times3(i)-pre_state_time):wake_times3(i)+post_state_time)-f1(i);
end