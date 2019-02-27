% test_FT_interTrialCoherence
% http://www.fieldtriptoolbox.org/faq/itc

if 0 % original from the link above (use ft_freqsimulation)
	cfg = [];
	cfg.numtrl	= 50;
	cfg.fsample     = 1000;
	cfg.trllen      = 1;
	
	cfg.s1.freq     = 30;
	cfg.s1.phase    =  'random'; % or 0
	cfg.s1.ampl     = 1;
		
% 	cfg.s2.freq = 80;
% 	cfg.s2.ampl = 1;
% 	cfg.s2.phase = 0;
% 	
 	cfg.noise.ampl = 1;

	cfg.method	= 'superimposed';
	cfg.output	= 'mixed';
	
	data = ft_freqsimulation(cfg); % simulate some data
	
else % create data "manually"
	fsample		= 1000;
	fsignal		= 30; % Hz
	nsamples	= 1000;
	taper		= 'hanning';
	ntrials		= 50;
	ampl1		= 1;
	snr		= 10; % inverse of noise amplitude
	phase_ran	= 2*pi/100; % 2*pi/1;  % max randomness: 2*pi/1

	phase1 = 2*pi * (zeros(1,ntrials) + phase_ran*rand(1,ntrials));

	
	data = [];
	data.label = {'ch1'};
	for i=1:ntrials
		data.time{i} = (1:nsamples)/fsample;
		data.trial{i} = [
			ampl1*cos(fsignal * 2*pi * data.time{i} + phase1(i)) + (1/snr) * randn(1,nsamples)/sqrt(2)
			];
	end
	
end

figure
cfg.channel = data.label{1};
d1 = ft_selectdata(cfg,data);

plot(data.time{1}, cell2mat(d1.trial')); title(sprintf('channel %s', data.label{1}));

cfg = [];

% cfg.method = 'mtmfft';
% cfg.output = 'pow';
% cfg.pad    = 'maxperlen';
% cfg.foilim = [0 100];
% cfg.taper  = 'hanning';

cfg.method = 'wavelet'; % wavelet
cfg.toi    = [0:0.01:1];
cfg.output = 'fourier';
cfg.foilim	= [0 100];
freq = ft_freqanalysis(cfg, data);

% figure
% switch cfg.output
% 	case 'fourier'
% 		plot(freq.freq, abs(squeeze(freq.fourierspctrm(1,1,:))), 'b-');
% 		sum(abs(freq.fourierspctrm(1,1,:)))
% 	case 'pow'
% 		% semilogy(freq.freq, freq.powspctrm(1,:), 'b-');
% 		plot(freq.freq, freq.powspctrm(1,:), 'b-');
% 		sum(freq.powspctrm(1,:)) % should be equal to 2*(A^2/4)
% end


% make a new FieldTrip-style data structure containing the ITC
% copy the descriptive fields over from the frequency decomposition

itc = [];
itc.label     = freq.label;
itc.freq      = freq.freq;
itc.time      = freq.time;
itc.dimord    = 'chan_freq_time';

F = freq.fourierspctrm;   % copy the Fourier spectrum
N = size(F,1);           % number of trials

% compute inter-trial phase coherence (itpc)
itc.itpc      = F./abs(F);         % divide by amplitude
itc.itpc      = sum(itc.itpc,1);   % sum angles
itc.itpc      = abs(itc.itpc)/N;   % take the absolute value and normalize
itc.itpc      = squeeze(itc.itpc); % remove the first singleton dimension

% compute inter-trial linear coherence (itlc)
itc.itlc      = sum(F) ./ (sqrt(N*sum(abs(F).^2)));
itc.itlc      = abs(itc.itlc);     % take the absolute value, i.e. ignore phase
itc.itlc      = squeeze(itc.itlc); % remove the first singleton dimension

% Finally we can plot it, just like a regular time-frequency representation

itc

figure
subplot(2, 1, 1);
imagesc(itc.time, itc.freq, squeeze(itc.itpc(:,:)));
% imagesc(itc.time, itc.freq, squeeze(itc.itpc(1,:,:))); 
axis xy
title('inter-trial phase coherence');
subplot(2, 1, 2);
imagesc(itc.time, itc.freq, squeeze(itc.itlc(:,:)));
% imagesc(itc.time, itc.freq, squeeze(itc.itlc(1,:,:)));
axis xy
title('inter-trial linear coherence');