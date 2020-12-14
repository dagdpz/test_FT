% function test_FT_spike_field_coherence
% http://www.fieldtriptoolbox.org/tutorial/spikefield
clear, clc, close; 

cfg = [];
cfg.numtrl	= 5;
cfg.fsample     = 1000; % Hz
cfg.trllen      = 1; % s

cfg.s1.freq     = 10; % Hz
cfg.s1.phase    = 0; %'random'; % 'random' or 0
cfg.s1.ampl     = 1;

cfg.s2.freq = 25; % Hz
cfg.s2.ampl = 0.5;
cfg.s2.phase = 'random';

cfg.s3.freq = 60; % Hz
cfg.s3.ampl = 0.2;
cfg.s3.phase = 0;

cfg.noise.ampl = 0.2;

cfg.method	= 'superimposed';
cfg.output	= 'all'; % mixed or all

data = ft_freqsimulation(cfg);

lfp_freq = [cfg.s1.freq cfg.s2.freq cfg.s3.freq];
lfp_amp =  [cfg.s1.ampl cfg.s2.ampl cfg.s3.ampl];

% figure('Name', 'Simulated LFP');
% plot(data.time{1}, data.trial{1}(1, 1:1000));
% title('Simulated LFP signal')

% spikes
spikeRate = [80]; % Hz, overall spike rate
spikePhaseFreq = 1; % 1 or 2 or 3, s1 or s2 or s3 components
spikePhaseMean = [-pi]; % align to trough
spikePhaseStd  = [pi/60]; % rad, spread around spikePhaseMean
spikeLockProb = 0.5; % probability of spike to be locked, the phase-locked firing rate would be max(spikeRate*spikeLockProb,lfp_freq(spikePhaseFreq))

% windos size 

lockedSpikes = zeros(cfg.numtrl,cfg.trllen*cfg.fsample);
Spikes = test_FT_simulate_spike_train(spikeRate*(1-spikeLockProb),cfg.trllen,cfg.numtrl); % non-locked spikes

extraSpikes = (spikeRate*spikeLockProb) - lfp_freq(spikePhaseFreq); % how many spikes will go into bursts

for t = 1:cfg.numtrl,
	% for each trial, find the phase (of the trough) of the corresponding freq component s
	 s = data.trial{t}(spikePhaseFreq+1,:);
	 taxis = data.time{t};
	 p = angle(hilbert(s));
	 trough_idx = find(diff(p)<-5); % find indices of the troughs
	 
     switch spikeRate*spikeLockProb <= lfp_freq(spikePhaseFreq);
         case 1 
             % If there is less locked spikes than cycles this line randomly
             % selects where to put them
             lockedSpikes_idx = trough_idx(rand(size(trough_idx)) <= spikeRate*spikeLockProb/lfp_freq(spikePhaseFreq));

             % add some variability to relative phase of spikes and troughs
             lockedSpikes_idx = fix(lockedSpikes_idx+randn(size(lockedSpikes_idx))*spikePhaseStd*(cfg.fsample/lfp_freq(spikePhaseFreq)));
             lockedSpikes_idx = ig_limit_range_min_max(lockedSpikes_idx,1,cfg.fsample*cfg.trllen);
             lockedSpikes(t,lockedSpikes_idx)=1;
             Spikes(t,lockedSpikes_idx) = 1; % add locked spikes to spikes
             data.trial{t}(6, :) = Spikes(t,:); 
         
         case 0
             % more locked spikes than cycles
             spikePerCycleAv = ceil(spikeRate*spikeLockProb / lfp_freq(spikePhaseFreq));
             % eto ya, kostyl`, but in fact i don`t know how to limit the
             % number of spikes per cycle
             if spikePerCycleAv > 4;
                 spikePerCycleAv = 4;
             end
             
             switch spikePhaseFreq
                 case 1
                     int = 5; % for each side
                 case 2
                     int = 3;
                 case 3
                     int = 2;
             end
             
             lockedSpikesBurstID = [];
             for trough = 1:numel(trough_idx);
                % generating few spikes in the interval around trough
                [low_r, up_r] = deal(trough_idx(trough) - int, trough_idx(trough) + int);
%                 low_r =  trough_idx(trough) - 5;
%                 up_r =  trough_idx(trough) + 5;
                range_id = round(low_r + (up_r - low_r).*rand(1,spikePerCycleAv));
                range_id = ig_limit_range_min_max(range_id, 1, cfg.fsample*cfg.trllen);
                lockedSpikesBurstID = [lockedSpikesBurstID, range_id];
             end
             lockedSpikes(t,lockedSpikesBurstID)=1;
             Spikes(t,lockedSpikesBurstID) = 1; % add locked spikes to spikes
             data.trial{t}(6, :) = Spikes(t,:); 
     end
end

disp(sprintf('spikeRate %d Hz, spikeLockProb %.2f, mean rate %.2f, mean locked rate %.2f %d spikes',spikeRate,spikeLockProb,mean(sum(Spikes,2)/cfg.trllen),mean(sum(lockedSpikes,2)/cfg.trllen),sum(sum(Spikes))));


data.label = {'lfp1', 's1', 's2', 's3', 'noise', 'spk1'};

MS_test_FT_plotTrial( 1, data, spikePhaseFreq, spikeRate, spikeLockProb, cfg.trllen, lockedSpikes, Spikes )
