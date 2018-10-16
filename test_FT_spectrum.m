% test_FT_spectrum
% run test_FT_initialize first, but only once

% http://www.fieldtriptoolbox.org/example/effects_of_tapering_for_power_estimates_in_the_frequency_domain
% http://www.gaussianwaves.com/2013/12/computation-of-power-of-a-signal-in-matlab-simulation-and-verification/

% 1. Spectrum
cfg = [];
cfg.method  = 'superimposed';
cfg.output = 'mixed';
cfg.fsample = 1000;
cfg.numtrl  = 1;
cfg.trllen  = 1; % s
cfg.s1.freq = 30;
cfg.s1.ampl = 3;
cfg.s1.phase = 0;
cfg.s2.freq = 80;
cfg.s2.ampl = 2; % use 0 to have only one signal
cfg.s2.phase = 0;
cfg.noise.ampl = 2;

data = ft_freqsimulation(cfg);
size(data)
figure
plot(data.time{1}, data.trial{1}(1,:))
title('Raw signal');
xlabel('Time (s)');
ylabel('Amplitude');

cfg        = [];
cfg.method = 'mtmfft'; % mtmfft | wavelet
cfg.output = 'fourier'; % fourier | pow
cfg.pad    = 'maxperlen';
cfg.foilim = [0 100];
cfg.taper  = 'hanning'; % hanning | dpss
cfg.taper     = 'dpss'; cfg.tapsmofrq = 1;   % i.e. no real smoothing

if strcmp(cfg.method,'wavelet'),
	cfg.toi    = 0.5;
	cfg.output = 'fourier';
end

freq       = ft_freqanalysis(cfg, data);
figure
switch cfg.output
	case 'fourier'
		plot(freq.freq, abs(squeeze(freq.fourierspctrm(1,1,:))), '.b-');
		sum(abs(freq.fourierspctrm(1,1,:))) 
	case 'pow'
		% semilogy(freq.freq, freq.powspctrm(1,:), 'b-');
		plot(freq.freq, freq.powspctrm(1,:), '.b-');
		sum(freq.powspctrm(1,:)) % should be equal to 2*(A^2/4)
end
title(sprintf('Spectrum, method %s output %s',cfg.method,cfg.output));
xlabel('Frequency (Hz)');
ylabel('Amplitude');

% 2. TFR
% http://www.fieldtriptoolbox.org/walkthrough

cfg             = [];
cfg.channel     = 'mix';
cfg.method      = 'mtmconvol'; % 'mtmconvol'; 'mtmfft'; 'wavelet';
cfg.taper       = 'dpss';  %'dpss', 'hanning';

cfg.foi         = [5:1:100]; % analysis ... to ... Hz in steps of ... Hz
cfg.output      = 'pow'; % 'pow'; % 'fourier';  % 'powandcsd'
cfg.keeptrials  = 'yes';
% cfg.padtype     ='edge';  % 'zero', 'mean', 'localmean', 'edge', 'mirror', 'nan' or 'remove'
% cfg.pad         = 1;
cfg.toi         = [0:0.025:1]; % s

switch cfg.method
	
	case 'mtmconvol'
		% cfg.t_ftimwin   = 5./cfg.foi; % CYCLES
		cfg.t_ftimwin   = ones(length(cfg.foi),1).* 0.250; % s, FIXED WIN
		
		switch cfg.taper
			
			case 'dpss'
				% cfg.tapsmofrq   = ones(length(cfg.foi),1).* 7; % FIXED
				% cfg.tapsmofrq	=  cfg.foi*0.3;
				cfg.tapsmofrq	=  [3*ones(length(find(cfg.foi<30)),1)' cfg.foi(cfg.foi>=30)*0.1];
			case 'hanning'
						
		end
		
	case 'wavelet'
		cfg.width      = 10;
end


TFR = ft_freqanalysis(cfg, data);

TFR_powspctrm= squeeze(TFR.powspctrm(1,1,:,:));
figure('Position',[100 100 1100 500]);
subplot(1,2,1)
% imagesc(TFR.time, TFR.freq, log10(TFR_powspctrm));
imagesc(TFR.time, TFR.freq, (TFR_powspctrm));
axis xy;
title('Time-frequency');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar

idx = find(TFR.time==0.5);
TFR_slice = squeeze(TFR_powspctrm(:,idx));
subplot(1,2,2)
plot(TFR.freq,TFR_slice);
title('Slice through 0.5 s');
xlabel('Frequency (Hz)');
ylabel('Amplitude');