% This script loops through different firing rate and plots ppc and coh 
% Based on test_FT_spike_field_coherence.m, check all the descriptions
% there
clear; close all; 

rates = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9]; %choose rates for looping
rates_labels = {'0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9'};
count = 1;

% parameters configuration for Chronux
params.Fs	=1000; % sampling frequency
params.fpass	=[1 100]; % band of frequencies to be kept
params.tapers	=[2 4]; % taper parameters
params.pad	=2; % pad factor for fft
params.err	=[2 0.05];
params.trialave =1;

for x=rates

% LFP signal is generated and kept constant
cfg = [];
cfg.numtrl	= 20;
cfg.fsample     = 1000; % Hz
cfg.trllen      = 1; % s

cfg.s1.freq     = 10; % Hz
cfg.s1.phase    = 'random'; % 'random' or 0
cfg.s1.ampl     = 1;

cfg.s2.freq = 25; % Hz
cfg.s2.ampl = 0.5;
cfg.s2.phase = 0;

cfg.s3.freq = 60; % Hz
cfg.s3.ampl = 0.2;
cfg.s3.phase = 0;

cfg.noise.ampl = 0.2;

cfg.method	= 'superimposed';
cfg.output	= 'all'; % mixed or all

data = ft_freqsimulation(cfg);
   
lfp_freq = [cfg.s1.freq cfg.s2.freq cfg.s3.freq];
lfp_amp =  [cfg.s1.ampl cfg.s2.ampl cfg.s3.ampl];

% spikes
spikeRate = [20]; % Hz, overall spike rate
spikePhaseFreq = 1; % 1 or 2 or 3, s1 or s2 or s3 components
spikePhaseMean = [-pi]; % align to trough
spikePhaseStd  = [pi/60]; % rad, spread around spikePhaseMean
spikeLockProb = x; % probability of spike to be locked, the phase-locked firing rate would be max(spikeRate*spikeLockProb,lfp_freq(spikePhaseFreq))

lockedSpikes = zeros(cfg.numtrl,cfg.trllen*cfg.fsample);
Spikes = test_FT_simulate_spike_train(spikeRate*(1-spikeLockProb),cfg.trllen,cfg.numtrl); % non-locked spikes
for t = 1:cfg.numtrl,
	% for each trial, find the phase (of the trough) of the corresponding freq component s
	 s = data.trial{t}(spikePhaseFreq+1,:);
	 taxis = data.time{t};
	 p = angle(hilbert(s));
	 trough_idx = find(diff(p)<-6); % find indices of the troughs
	 
	 lockedSpikes_idx = trough_idx(rand(size(trough_idx)) <= spikeRate*spikeLockProb/lfp_freq(spikePhaseFreq));
	 
	 % add some variability to relative phase of spikes and troughs
	 lockedSpikes_idx = fix(lockedSpikes_idx+randn(size(lockedSpikes_idx))*spikePhaseStd*(cfg.fsample/lfp_freq(spikePhaseFreq)));
	 lockedSpikes_idx = ig_limit_range_min_max(lockedSpikes_idx,1,cfg.fsample*cfg.trllen);
	 lockedSpikes(t,lockedSpikes_idx)=1;
	 Spikes(t,lockedSpikes_idx) = 1; % add locked spikes to spikes
	 data.trial{t} = [data.trial{t}; Spikes(t,:)]; 
	 
	 if 0 % plot lfp and spikes (for illustration and debugging)
		 spikes_idx = find(Spikes(t,:));
		 rand_rgb = rand(1,3);
		 plot(taxis,s,'Color',rand_rgb); hold on;
		 plot(taxis(spikes_idx),-1*lfp_amp(spikePhaseFreq),'rx','MarkerEdgeColor',rand_rgb); % all spikes
		 if ~isempty(lockedSpikes_idx)
			plot(taxis(lockedSpikes_idx),-1*lfp_amp(spikePhaseFreq),'ro','MarkerEdgeColor',rand_rgb); % locked spikes
		 end
		 pause;
	 end
	
end

disp(sprintf('spikeRate %d Hz, spikeLockProb %.2f, mean rate %.2f, mean locked rate %.2f %d spikes',spikeRate,spikeLockProb,mean(sum(Spikes,2)/cfg.trllen),mean(sum(lockedSpikes,2)/cfg.trllen),sum(sum(Spikes))));

data.label = {'lfp1', 's1', 's2', 's3', 'noise', 'spk1'};
% data_all = ft_appendspike([],data, Spikes)

disp('Spike-triggered spectrum, mtmconvol');
% method 2
cfg           = [];
cfg.method    = 'mtmconvol';
cfg.foi       = 5:1:100;
cfg.t_ftimwin = 5./cfg.foi; % 5 cycles per frequency
cfg.taper     = 'hanning';
cfg.spikechannel = data.label{6};
cfg.channel      = data.label{1};
stsConvol     = ft_spiketriggeredspectrum(cfg, data);
 
% compute the statistics on the phases
cfg               = [];
cfg.method        = 'ppc0'; % compute the Pairwise Phase Consistency, can also be plv, ang, ppc1, ppc2, ral
cfg.spikechannel  = stsConvol.label{1};
cfg.channel       = stsConvol.lfplabel(1); % selected LFP channels
cfg.avgoverchan   = 'unweighted'; % weight spike-LFP phases irrespective of LFP power
cfg.timwin        = 'all'; % compute over all available spikes in the window
cfg.latency       = 'maxperiod'; % [0 nanmax(stsConvol.trialtime(:))]; %
statSts           = ft_spiketriggeredspectrum_stat(cfg,stsConvol);

ppc.freq = statSts.freq
ppc.coef(count, :) = statSts.ppc0

%% Compute Coherence (coh) through ft_freqanalysis and ft_connectivityanalysis

% for convenience select only LFP and spike data
cfg = []
cfg.channel = {'lfp1', 'spk1'}
data_selected = ft_selectdata(cfg, data)

cfg = [];
cfg.method	    = 'mtmfft';
cfg.taper	    = 'hanning'; % 'hanning' | 'dpss'
cfg.tapsmofrq	= 4; % for dpss
cfg.output      = 'powandcsd'; % powandcsd | fourier (for amplcorr or powcorr)
cfg.foilim      = [0 100];
freq            = ft_freqanalysis(cfg, data_selected);

cfg = [];
cfg.method	= 'coh'; % coh plv ppc amplcorr powcorr  
cfg.complex	= 'complex'; % abs complex

conn = ft_connectivityanalysis(cfg, freq);

coh.freq = conn.freq;
coh.coef(count, :) = abs(conn.cohspctrm);

% Calculating coherency with Chronux
% LFP: time x trials
cfg			= [];
cfg.trials		= 'all';
cfg.channel		= data.label{1};
data_lfp		= ft_selectdata(cfg, data);
chr_data1		= test_FT_fieldtrip2chronux(data_lfp,'lfp');

cfg.channel		= data.label{6};
data_spikes		= ft_selectdata(cfg, data);
chr_data2		= test_FT_fieldtrip2chronux(data_spikes,'spikes');

[C,phi,S12,S1,S2,f,zerosp,confC,phistd,Cerr]=coherencycpt(chr_data1,chr_data2,params,0);

chronux_coh.freq = f;
chronux_coh.coef(count, :) = C;

count = count + 1;
end

%% 
figure('Name', 'Coh and PPC comparison with different FR') % add LFP simulation parameters here
title('whatever')
subplot(1,3,1)
plot(coh.freq, coh.coef)
xlabel('Frequency')
ylabel('Coh')
title('Coherence using FT')

subplot(1,3,2)
plot(chronux_coh.freq, chronux_coh.coef)
xlabel('Frequency')
ylabel('Coh')
title('Coh using Chronux')

subplot(1,3,3)
plot(ppc.freq, ppc.coef)
xlabel('Frequency')
ylabel('Ppc0')
title('PPC using FT')
legend(rates_labels)

%% Compare coh in FT and Cronux
figure('Name', 'FT vs Chronux coh comparison');
freq_man = linspace(1.2207,100,length(chronux_coh.coef(5, :))) %doublecheck how it is constructed and whether introduces shift
for m = 1:length(rates)
    subplot(3,3,m)
    plot(coh.coef(m, :), 'r', 'DisplayName', 'FT'); hold on; plot(freq_man, chronux_coh.coef(m, :))
    xlim([0 100])
    ylim([0 1])
    legend('FT', 'Chronux')
    title(['Spiking probability is ', rates_labels(m)])
end

% legend, x up to 100, title, y up to 1

%% Plot the output of Chronux - what are the other two? 

% figure('Name','Chronux coherencycpt'); 
% subplot(311); plot(f,C); hold on; plot(f,Cerr(1,:),'r:'); plot(f,Cerr(2,:),'r:');
% subplot(312); plot(f,10*log10(S1));
% subplot(313); plot(f,10*log10(S2));
% xlabel('frequency')