function [ output_args ] = plotTrial( which_trial, data, spikePhaseFreq, spikeRate, spikeLockProb, trllen, lockedSpikes, Spikes )
% Use this function to take a look at how simulated data looks like

subplot(4, 1, 1);
plot(data.time{1}, data.trial{which_trial}(1, :)); % mixed signal
title('Mixed LFP signal');

subplot(4,1,2);
plot(data.time{1}, data.trial{which_trial}(spikePhaseFreq+1, :)); % LFP channel with locked spikes
title('LFP component with locked spikes');

subplot(4,1,3);
plotRaster(logical(lockedSpikes), data.time{1});
title(sprintf('Locked spikes; FR %d Hz, spikeLockProb %.2f, mean locked rate %.2f, %d locked spikes', ... 
    spikeRate,spikeLockProb,mean(sum(lockedSpikes,2)/trllen),sum(sum(lockedSpikes))));

subplot(4,1,4);
plotRaster(Spikes, data.time{1});
title(sprintf('All spikes; FR %d Hz, spikeLockProb %.2f, mean rate %.2f, %d spikes',...
    spikeRate,spikeLockProb,mean(sum(Spikes,2)/trllen),sum(sum(Spikes))));

end

