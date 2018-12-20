% This script loops through different firing rate and plots ppc and coh 
% Based on test_FT_spike_field_coherence.m, check all the descriptions
% there
clear; close all; clc; 

%% what to run in this script
compute_coh_chronux = 0; % 0 for no
compute_coh_ft = 0;
compute_ppc0 = 0;
plot_separately = 0;

comparePPC = 1;
plot_compared_ppc = 1;

%% configs for LFP mixed signal generation with FieldTrip
cfg = [];
cfg.numtrl	= 20;
cfg.fsample     = 1000; % Hz
cfg.trllen      = 1; % s

cfg.s1.freq     = 10; % Hz
cfg.s1.phase    = 0; % 'random' or 0
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

%% 
spikeLockProbs = 0:0.1:1; %choose rates for looping
spikeLockProbs_labels = {'0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1'};

spikeRates = 10:10:100;

ppc0 = zeros(length(spikeLockProbs), length(spikeRates));
ppc1 = zeros(length(spikeLockProbs), length(spikeRates));
coh_ft = zeros(length(spikeLockProbs), length(spikeRates));
coh_chronux = zeros(length(spikeLockProbs), length(spikeRates));

%%
for lockv=1:length(spikeLockProbs);

% LFP signal is generated and kept constant
data = ft_freqsimulation(cfg);
   
lfp_freq = [cfg.s1.freq cfg.s2.freq cfg.s3.freq];
lfp_amp =  [cfg.s1.ampl cfg.s2.ampl cfg.s3.ampl];

spikePhaseFreq = 1; % 1 or 2 or 3, s1 or s2 or s3 components !!! IMPORTANT - CHANGE FOR LOCKING FREQ - 1 for 10hz, 2 for 25, 3 for 60
spikePhaseMean = [-pi]; % align to trough
spikePhaseStd  = [pi/60]; % rad, spread around spikePhaseMean
spikeLockProb = spikeLockProbs(lockv); % probability of spike to be locked, the phase-locked firing rate would be max(spikeRate*spikeLockProb,lfp_freq(spikePhaseFreq))

    for ratev = 1:length(spikeRates);
    % spikes
    spikeRate = [spikeRates(ratev)]; % Hz, overall spike rate


    lockedSpikes = zeros(cfg.numtrl,cfg.trllen*cfg.fsample);
    Spikes = test_FT_simulate_spike_train(spikeRate*(1-spikeLockProb),cfg.trllen,cfg.numtrl); % non-locked spikes
    for t = 1:cfg.numtrl,
        % for each trial, find the phase (of the trough) of the corresponding freq component s
         s = data.trial{t}(spikePhaseFreq+1,:);
         taxis = data.time{t};
         p = angle(hilbert(s));
         trough_idx = find(diff(p)<-5); % find indices of the troughs
	 
         lockedSpikes_idx = trough_idx(rand(size(trough_idx)) <= spikeRate*spikeLockProb/lfp_freq(spikePhaseFreq));
	 
         % add some variability to relative phase of spikes and troughs
         lockedSpikes_idx = fix(lockedSpikes_idx+randn(size(lockedSpikes_idx))*spikePhaseStd*(cfg.fsample/lfp_freq(spikePhaseFreq)));
         lockedSpikes_idx = ig_limit_range_min_max(lockedSpikes_idx,1,cfg.fsample*cfg.trllen);
         lockedSpikes(t,lockedSpikes_idx)=1;
         Spikes(t,lockedSpikes_idx) = 1; % add locked spikes to spikes
         data.trial{t}(6, :) = Spikes(t,:); 
	 
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

    %% Compute ppc with mtmconvol
    % method 2
    if compute_ppc0 == 1
        cfg_ppc           = [];
        cfg_ppc.method    = 'mtmconvol';
        cfg_ppc.foi       = lfp_freq(spikePhaseFreq);
        cfg_ppc.t_ftimwin = 5./cfg_ppc.foi; % 5 cycles per frequency
        cfg_ppc.taper     = 'hanning';
        cfg_ppc.spikechannel = data.label{6};
        cfg_ppc.channel      = data.label{1};
        stsConvol     = ft_spiketriggeredspectrum(cfg_ppc, data);

        % compute the statistics on the phases
        cfg_ppc               = [];
        cfg_ppc.method        = 'ppc0'; % compute the Pairwise Phase Consistency, can also be plv, ang, ppc1, ppc2, ral
        cfg_ppc.spikechannel  = stsConvol.label{1};
        cfg_ppc.channel       = stsConvol.lfplabel(1); % selected LFP channels
        cfg_ppc.avgoverchan   = 'unweighted'; % weight spike-LFP phases irrespective of LFP power
        cfg_ppc.timwin        = 'all'; % compute over all available spikes in the window
        cfg_ppc.latency       = 'maxperiod'; % [0 nanmax(stsConvol.trialtime(:))]; %
        statSts           = ft_spiketriggeredspectrum_stat(cfg_ppc,stsConvol);

        ppc0(lockv, ratev) = statSts.ppc0;
        
    end
    
    %% Compute Coherence (coh) through ft_freqanalysis and ft_connectivityanalysis
    % for convenience select only LFP and spike data
    if compute_coh_ft == 1
        cfg_coh = [];
        cfg_coh.channel = {'lfp1', 'spk1'};
        data_selected = ft_selectdata(cfg_coh, data);

        cfg_coh = [];
        cfg_coh.method	    = 'mtmfft';
        cfg_coh.taper	    = 'hanning'; % 'hanning' | 'dpss'
        cfg_coh.tapsmofrq	= 4; % for dpss
        cfg_coh.output      = 'powandcsd'; % powandcsd | fourier (for amplcorr or powcorr)
        cfg_coh.foilim      = [0 100];
        freq            = ft_freqanalysis(cfg_coh, data_selected);

        cfg_coh = [];
        cfg_coh.method	= 'coh'; % coh plv ppc amplcorr powcorr  
        cfg_coh.complex	= 'complex'; % abs complex

        conn = ft_connectivityanalysis(cfg_coh, freq);
        ind_freq = find(conn.freq == lfp_freq(spikePhaseFreq));

        coh_ft(lockv, ratev) = abs(conn.cohspctrm(ind_freq));
    end
    
    %% Compute coh using chronux
    if compute_coh_chronux == 1;
        params.Fs	=1000; % sampling frequency
        params.fpass	=[1 lfp_freq(spikePhaseFreq)]; % band of frequencies to be kept
        params.tapers	=[2 4]; % taper parameters
        params.pad	=2; % pad factor for fft
        params.err	=[2 0.05];
        params.trialave =1; 

    	% Calculating coherency with Chronux
        % LFP: time x trials
        cfg_ch			= [];
        cfg_ch.trials		= 'all';
        cfg_ch.channel		= data.label{1};
        data_lfp		= ft_selectdata(cfg_ch, data);
        chr_data1		= test_FT_fieldtrip2chronux(data_lfp,'lfp');

        cfg_ch.channel		= data.label{6};
        data_spikes		= ft_selectdata(cfg_ch, data);
        chr_data2		= test_FT_fieldtrip2chronux(data_spikes,'spikes');

        [C,phi,S12,S1,S2,f,zerosp,confC,phistd,Cerr]=coherencycpt(chr_data1,chr_data2,params,0);

        coh_chronux(lockv, ratev) = C(end);
    end
    
    %% Comparing ppc0 and ppc1 from ft_spiketriggeredspectrum_stat
    
    if comparePPC == 1;
        cfg_ppc           = [];
        cfg_ppc.method    = 'mtmconvol';
        cfg_ppc.foi       = lfp_freq(spikePhaseFreq);
        cfg_ppc.t_ftimwin = 5./cfg_ppc.foi; % 5 cycles per frequency
        cfg_ppc.taper     = 'hanning';
        cfg_ppc.spikechannel = data.label{6};
        cfg_ppc.channel      = data.label{1};
        stsConvol     = ft_spiketriggeredspectrum(cfg_ppc, data);

        % compute the statistics on the phases
        cfg_ppc               = [];
        cfg_ppc.spikechannel  = stsConvol.label{1};
        cfg_ppc.channel       = stsConvol.lfplabel(1); % selected LFP channels
        cfg_ppc.avgoverchan   = 'unweighted'; % weight spike-LFP phases irrespective of LFP power
        cfg_ppc.timwin        = 'all'; % compute over all available spikes in the window
        cfg_ppc.latency       = 'maxperiod'; % [0 nanmax(stsConvol.trialtime(:))]; %
        
        % ppc0
        cfg_ppc.method        = 'ppc0'; % compute the Pairwise Phase Consistency, can also be plv, ang, ppc1, ppc2, ral
        statStsppc0           = ft_spiketriggeredspectrum_stat(cfg_ppc,stsConvol);
        ppc0(lockv, ratev)    = statStsppc0.ppc0;
        
        %ppc1
        cfg_ppc.method        = 'ppc1'; % compute the Pairwise Phase Consistency, can also be plv, ang, ppc1, ppc2, ral
        statStsppc1           = ft_spiketriggeredspectrum_stat(cfg_ppc,stsConvol);
        ppc1(lockv, ratev) = statStsppc1.ppc1;
    end    
    
    end
end


%% Now output for certain frequency is in ppc, coh_ft, and coh_chronux
%  Rows are different locking probabilities, columns are different firing
%  rates

keys.path_to_save = sprintf('%s\\plots', pwd);

if exist(sprintf('%s\\plots', pwd)) == 0
    mkdir(sprintf('%s\\plots', pwd))
end

%% Plotting ppc0, coh_ft, coh_chronux separately
if plot_separately == 1;
    coef.ppc0 = ppc0;
    coef.coh_ft = coh_ft;
    coef.coh_chronux = coh_chronux;

    names = fieldnames(coef);
    rand_rgb = rand(numel(spikeLockProbs),3);

    for c = 1:numel(names);
        coef_tmp = getfield(coef, names{c});

        FigName_short1= sprintf('%s_%dHz_LFP', names{c}, lfp_freq(spikePhaseFreq)); %name of the file
        FigName_Long1=sprintf('%s value with locking to %d Hz LFP', names{c}, lfp_freq(spikePhaseFreq)); % overall title of the figure
        set(0,'DefaultTextInterpreter','none');
        wanted_papersize=[40 25];
        h1=figure('outerposition',[0 0 1900 1200],'name', FigName_Long1);
        figure(h1)

        for el = 1:numel(spikeLockProbs);
            figure(c); hold on
            plot(spikeRates, coef_tmp(el, :), 'o-', 'Color', rand_rgb(el,:), 'MarkerFaceColor', rand_rgb(el,:))
            xlabel('spike rate (Hz)'); ylabel(sprintf('%s', names{c})); ylim([0 1])
            lgd1 = legend('0.0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1.0');
            title(lgd1, 'spikeLockProb');
            %title('Pairwise phase consistency (ppc) at 10 Hz')
            grid on
        end    

        set(figure(h1), 'Paperunits','centimeters','PaperSize', wanted_papersize,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_papersize])
        mtit(figure(h1),  [FigName_Long1 ], 'xoff', 0, 'yoff', 0.05, 'color', [0 0 0], 'fontsize', 14,'Interpreter', 'none');
        export_fig(figure(h1), [keys.path_to_save filesep FigName_short1], '-pdf','-transparent') % pdf by run
    end
end

%% Plotting ppc0 and ppc1 to compare
if plot_compared_ppc == 1;
    rand_rgb = rand(numel(spikeLockProbs),3);
    
    FigName_short1= sprintf('PPC_comparison_%dHz_LFP', lfp_freq(spikePhaseFreq)); %name of the file
    FigName_Long1=sprintf('PPC value with locking to %d Hz LFP', lfp_freq(spikePhaseFreq)); % overall title of the figure
    set(0,'DefaultTextInterpreter','none');
    wanted_papersize=[40 25];
    h1=figure('outerposition',[0 0 1900 1200],'name', FigName_Long1);
    figure(h1)
    
    for el = 1:numel(spikeLockProbs);
        figure(1); hold on
        plot(spikeRates, ppc0(el, :), 'o-', 'Color', rand_rgb(el,:), 'MarkerFaceColor', rand_rgb(el,:), 'DisplayName', 'ppc0');
        plot(spikeRates, ppc1(el, :), 'o--', 'Color', rand_rgb(el,:), 'MarkerFaceColor', rand_rgb(el,:), 'DisplayName', 'ppc1')
        xlabel('spike rate (Hz)'); ylabel('ppc'); ylim([0 1])
        lgd1 = legend('0.0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1.0');
        title(lgd1, 'spikeLockProb');
        grid on
    end
    
    set(figure(h1), 'Paperunits','centimeters','PaperSize', wanted_papersize,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_papersize])
    mtit(figure(h1),  [FigName_Long1 ], 'xoff', 0, 'yoff', 0.05, 'color', [0 0 0], 'fontsize', 14,'Interpreter', 'none');
    export_fig(figure(h1), [keys.path_to_save filesep FigName_short1], '-pdf','-transparent') % pdf by run
end