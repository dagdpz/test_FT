% test_FT_interTrial_variability
% Simulate EEG/LFP variability testing

% Useful info
% http://www.fieldtriptoolbox.org/_media/workshop/marseille_frequency.pdf
% see also http://www.fieldtriptoolbox.org/faq/itc

% 0. Need to run this once if fieldtrip is not already on the path
% addpath('D:\Sources\fieldtrip-20151228');
% ft_defaults


% 1. Manually create data, baseline (prestim) and response (poststim) periods, each 1 s
% (also can use ft_freqsimulation, in principle, but it is not as flexible)

fsample		= 1000;
nsamples	= 1000; % for each period (baseline/response)
taper		= 'hanning';
ntrials		= 20;

% baseline
b.freq1		= 10 + 0.0*randn(1,ntrials); % Hz, the bigger is multiplier, the larger is trial-to-trial variability, around base freq 
b.ampl1		= 2;
b.ampl1_ran	= 0.01*b.ampl1*randn(1,ntrials); % the bigger is multiplier, the larger is trial-to-trial variability 
b.freq2		= 30 + 0.0*randn(1,ntrials); % Hz
b.ampl2		= 1;
b.ampl2_ran	= 0.01*b.ampl2*randn(1,ntrials);

b.phase1_ran	= 2*pi/100; % [from 0 to max randomness: 2*pi/1]
b.phase1	= 2*pi * (zeros(1,ntrials) + b.phase1_ran*rand(1,ntrials));
b.phase2_ran	= 0; % [from 0 to max randomness: 2*pi/1]
b.phase2	= 2*pi * (zeros(1,ntrials) + b.phase2_ran*rand(1,ntrials));

% response
r.freq1		= 10 + 0.9*randn(1,ntrials); % Hz
r.ampl1		= 7;
r.ampl1_ran	= 0.25*r.ampl1*randn(1,ntrials); % the bigger is multiplier, the larger is trial-to-trial variability  
r.freq2		= 30 + 2*randn(1,ntrials); % Hz
r.ampl2		= 3;
r.ampl2_ran	= 0.25*r.ampl2*randn(1,ntrials); % the bigger is multiplier, the larger is trial-to-trial variability


r.phase1_ran	= 2*pi/10; % [from 0 to max randomness: 2*pi/1]
r.phase1	= 2*pi * (zeros(1,ntrials) + r.phase1_ran*rand(1,ntrials));
r.phase2_ran	= 0; % [from 0 to max randomness: 2*pi/1]
r.phase2	= 2*pi * (zeros(1,ntrials) + r.phase2_ran*rand(1,ntrials));

% IMPORTANT NOTE about phase: if response phase is selected independently from baseline, this would result in phase "reset" even if response phase is random
% therefore use baseline phase if no phase reset need to be modelled (uncomment below)
% r.phase1 = b.phase1; r.phase2 = b.phase2; % uncomment this for smooth phase from baseline to response

noiseampl	= 0.1; % additive Gaussian noise amplitude


data = [];
data.label = {'ch1'};
for i=1:ntrials
	data.time{i} = (-nsamples:nsamples-1)/fsample;
	data.trial{i} = [
		max(0,(b.ampl1+b.ampl1_ran(i)))*cos(b.freq1(i) * 2*pi * data.time{i}(1:nsamples) + b.phase1(i)) + ...
		max(0,(b.ampl2+b.ampl2_ran(i)))*cos(b.freq2(i) * 2*pi * data.time{i}(1:nsamples) + b.phase2(i)) + ...
		+ noiseampl * randn(1,nsamples)...
		max(0,(r.ampl1+r.ampl1_ran(i)))*cos(r.freq1(i) * 2*pi * data.time{i}(nsamples+1:2*nsamples) + r.phase1(i)) + ...
		max(0,(r.ampl2+r.ampl2_ran(i)))*cos(r.freq2(i) * 2*pi * data.time{i}(nsamples+1:2*nsamples) + r.phase2(i)) + ...
		+ noiseampl * randn(1,nsamples)
		];
end

% 2. Plot time domain data and spectrum
figure('Position',[100 100 800 800]);
subplot(2,1,1);
cfg.latency = 'prestim';
d_b = ft_selectdata(cfg,data);
plot(d_b.time{1}, cell2mat(d_b.trial'),'LineWidth',0.5); hold on

cfg.latency = 'poststim';
d_r = ft_selectdata(cfg,data);
plot(d_r.time{1}, cell2mat(d_r.trial'),'LineWidth',1); hold on


plot(data.time{1},mean(cell2mat(data.trial'),1),'k-','LineWidth',2); % average (evoked) response

title(sprintf('channel %s', data.label{1}));
xlabel('Time from stim onset (s)');

cfg = [];

cfg.method = 'mtmfft';
cfg.output = 'pow';
cfg.pad    = 'maxperlen';
cfg.foilim = [0 100];
cfg.taper  = 'hanning';

% OR
% cfg.method = 'wavelet'; 
% cfg.output = 'fourier';
% cfg.foilim = [0 100];
 
cfg.toi    = [-0.5]; % for wavelet
freq_b = ft_freqanalysis(cfg, d_b);
cfg.toi    = [0.5]; % for wavelet
freq_r = ft_freqanalysis(cfg, d_r);

subplot(2,1,2)
switch cfg.output
	case 'fourier'
		plot(freq_b.freq, abs(squeeze(freq_b.fourierspctrm(1,1,:))), 'b-'); hold on
		plot(freq_r.freq, abs(squeeze(freq_r.fourierspctrm(1,1,:))), 'r-');
	case 'pow'
		% or use semilogy
		plot(freq_b.freq, freq_b.powspctrm(1,:), 'b-'); hold on
		plot(freq_r.freq, freq_r.powspctrm(1,:), 'r-'); 
		legend({'prestim' 'poststim'});
end
xlabel('Hz');
legend({'prestim' 'poststim'});

% 3. Time-Frequency analysis (on combined [prestim poststim])
% http://www.fieldtriptoolbox.org/walkthrough

cfg             = [];
cfg.method      = 'mtmconvol'; % 'mtmconvol'; 'mtmfft'; 'wavelet';
cfg.taper       = 'dpss';  %'dpss', 'hanning';

cfg.foi         = [5:1:100]; % analysis ... to ... Hz in steps of ... Hz
cfg.output      = 'pow'; % 'pow'; % 'fourier';  % 'powandcsd'
cfg.keeptrials  = 'yes'; % save analysis of individual trials
% cfg.padtype     ='edge';  % 'zero', 'mean', 'localmean', 'edge', 'mirror', 'nan' or 'remove'
% cfg.pad         = 1;
cfg.toi         = [-1:0.025:1]; % s

switch cfg.method
	
	case 'mtmconvol'
		% cfg.t_ftimwin   = 5./cfg.foi; % CYCLES
		cfg.t_ftimwin   = ones(length(cfg.foi),1).* 0.250; % s, FIXED WIN
		
		switch cfg.taper
			
			case 'dpss' % good for high frequencies
				% cfg.tapsmofrq = ones(length(cfg.foi),1)*7; % FIXED
				% cfg.tapsmofrq = 5;
				% cfg.tapsmofrq	=  cfg.foi*0.3;
				cfg.tapsmofrq	=  [4*ones(length(find(cfg.foi<40)),1)' cfg.foi(cfg.foi>=40)*0.1]; % adaptive number of tapers, increasing with freq
				
			case 'hanning' % good for low frequencies
						
		end
		
	case 'wavelet'
		cfg.width      = 5;
end

TFR = ft_freqanalysis(cfg, data);

TFR_powspctrm_average = squeeze(mean(TFR.powspctrm(:,1,:,:),1)); % mean across all trials
figure('Position',[100 100 1300 600]);
subplot(1,2,1)
% imagesc(TFR.time, TFR.freq, log10(TFR_powspctrm)); 
% imagesc(TFR.time, TFR.freq, zscore(TFR_powspctrm_average,2)); % make sure to use zscore from NaN toolbox
imagesc(TFR.time, TFR.freq, TFR_powspctrm_average);
axis xy;
title('Time-frequency');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar

% Band-limited power - derived from above Time-Freq analysis
band(1,:)	= [5 15]; % Hz
band(2,:)	= [25 35];% Hz

analWin(1,:)	= [-0.8 -0.2];	% s, baseline (prestim)
analWin(2,:)	= [0.2 0.8];	% s, response (poststim)

NORMALIZE_BY_ZSCORE = 1; % 0 or 1, if 1 normalize each band by zscore across all concatinated trials, to disregard 1/f power law

% Pre-allocate
R(ntrials,size(band,1),size(analWin,1))=NaN; % Response measure [trials x band x analysis window]

subplot(1,2,2)
cmap = jet(size(band,1));

for b = 1:size(band,1), % for each band
	blp = squeeze( nanmean(TFR.powspctrm( :, 1, TFR.freq >= band(b,1) & TFR.freq <= band(b,2) ,:),3)); % Band-limited power
	
	if NORMALIZE_BY_ZSCORE
		conc_blp = reshape(blp,1,ntrials*size(blp,2));
		blp = (blp - mean( conc_blp ))/nanstd( conc_blp );
		zscore_string = 'zscored';
	else
		zscore_string = '';
	end
	
	hp(b,:) = plot(TFR.time, blp ,'Color',cmap(b,:)); hold on % plot all trials
	for aw = 1:size(analWin,1),
		R(:,b,aw) = mean(blp(:, TFR.time >= analWin(aw,1) & TFR.time <= analWin(aw,2)),2);
	end
	
end
legend(hp(:,1)',num2str(band(:,:)),'Location','best');
title(sprintf('Band-limited power %s', zscore_string));
xlabel('Time (s)');
ylabel('Power');

% Response variability analysis
figure('Position',[100 100 400 400]);
for b = 1:size(band,1), % for each band
	for aw = 1:size(analWin,1),
		plot(aw,mean(abs(R(:,b,aw) - mean(R(:,b,aw),1))),'o','MarkerSize',10,'Color',cmap(b,:)); hold on
	end	
end
xlabel('Analysis window');
title(['mean(abs(Power_{i} - mean(Power_{i})) ' zscore_string]);
set(gca,'Xlim',[0 size(analWin,1)+1],'Xtick',[1:1:size(analWin,1)]);

