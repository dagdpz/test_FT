% The Coherence spectrum
% http://www.fieldtriptoolbox.org/tutorial/fourier
% http://www.fieldtriptoolbox.org/example/coherence_snr

fsample		= 1000;
fsignal		= 10; % Hz
nsamples	= 1000;
taper		= 'hanning';
ntrials		= 50;
snr		= 100;
phase_ran	= 2*pi; % 2*pi/1;  % max randomness: 2*pi/1
ampl1_ran	= 0; % ampl variability ch1
ampl2_ran	= 0; % (independent) ampl variability ch2
ampl12_ran	= 0; % strength of correlation between two channels

% introduce a random phase difference between the two channels on each trial
phase1 = 2*pi * (zeros(1,ntrials)		     );
phase2 = 2*pi * (zeros(1,ntrials) + phase_ran*rand(1,ntrials));

% introduce a amplitude modulation between the two channels on each trial
ampl1  = 1 + ampl1_ran*randn(1,ntrials);
ampl2  = 1 + ampl2_ran*randn(1,ntrials);
ampl2  = ampl2 + ampl12_ran*ampl1;


data = [];
data.label = {'ch1', 'ch2'};
for i=1:ntrials
  data.time{i} = (1:nsamples)/fsample;
  data.trial{i} = [
    ampl1(i)*cos(fsignal * 2*pi * data.time{i} + phase1(i)) + (1/snr) * randn(1,nsamples)/sqrt(2);
    ampl2(i)*cos(fsignal * 2*pi * data.time{i} + phase2(i)) + (1/snr) * randn(1,nsamples)/sqrt(2);
    ];
end

figure
d1 = ft_selectdata(data,'channel','ch1');
d2 = ft_selectdata(data,'channel','ch2');

subplot(2,1,1); plot(data.time{1}, cell2mat(d1.trial')); title(sprintf('channel %s', data.label{1}));
subplot(2,1,2); plot(data.time{1}, cell2mat(d2.trial')); title(sprintf('channel %s', data.label{2}));

cfg = [];
cfg.method	= 'mtmfft';
cfg.taper	= 'hanning'; % 'hanning' | 'dpss'
cfg.tapsmofrq	= 4; % for dpss
cfg.output	= 'powandcsd'; % powandcsd | fourier (for amplcorr or powcorr)
cfg.foilim	= [0 100];
freq = ft_freqanalysis(cfg, data);


cfg = [];
cfg.method	= 'coh'; % coh plv ppc amplcorr powcorr  
cfg.complex	= 'complex'; % abs complex

% cfg.jackknife = 'yes';

conn = ft_connectivityanalysis(cfg, freq);

figure
hold on

switch cfg.method

	case 'coh'
		plot(conn.freq, abs(conn.cohspctrm));
	case 'plv'
		plot(conn.freq, abs(conn.plvspctrm));
	case 'ppc'
		plot(conn.freq, abs(conn.ppcspctrm));
	case 'amplcorr' % needs fourier
		plot(conn.freq, abs(squeeze(conn.amplcorrspctrm(1,2,:))));
	case 'powcorr' % needs fourier
		plot(conn.freq, abs(squeeze(conn.powcorrspctrm(1,2,:))));
end