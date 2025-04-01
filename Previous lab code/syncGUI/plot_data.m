
for ix = 1
figure
subplot(2,1,1)
plot(trials(ix).ts,trials(ix).frameReadOut)
title('Hamamatsu')
ylim([0 3])
xlim([0 5])
subplot(2,1,2)
plot(trials(ix).ts,trials(ix).baslerExposure)
title('Basler')
ylim([0 5])
xlim([0 5])
end

%%
for ix = 1:length(trials)
figure
subplot(3,1,1)
plot(trials(ix).ts,trials(ix).frameReadOut)
title('Hamamatsu')
ylim([0 3])
subplot(3,1,2)
plot(trials(ix).ts,trials(ix).stimulusTrigger)
title('Stimulus')
ylim([0 5])
subplot(3,1,3)
plot(trials(ix).ts,trials(ix).baslerExposure)
title('Basler')
ylim([0 5])
end