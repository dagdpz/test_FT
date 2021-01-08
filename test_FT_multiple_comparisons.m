% function test_FT_multiple_comparisons
% run test_FT_initialize first, but only once
% https://www.fieldtriptoolbox.org/workshop/oslo2019/introduction/
% https://www.fieldtriptoolbox.org/workshop/oslo2019/statistics/

numtrl = 100; % for each of two conditions
noise1 = 1;
noise2 = 1;
ampl_diff = 1;
alpha = 0.05;
n_samples = 200;

idx_diff = [160:170 180:181 185:195]; % indices when data1 ~= data2
idx_same = setdiff(1:n_samples,idx_diff);

% signals
cfg = [];
cfg.method  = 'superimposed';
cfg.output = 'mixed';
cfg.fsample = n_samples;
cfg.numtrl  = numtrl;
cfg.trllen  = 1; % s
cfg.s1.freq = 1;
cfg.s1.ampl = 5;
cfg.s1.phase = 0;
cfg.noise.ampl = noise1;

data1 = ft_freqsimulation(cfg);

cfg.s1.freq = 1;
cfg.s1.ampl = 5-ampl_diff;
cfg.s1.phase = 0;
cfg.noise.ampl = noise2;

data2 = ft_freqsimulation(cfg);

cfg.s1.freq = 1;
cfg.s1.ampl = 5;
cfg.s1.phase = 0;
cfg.noise.ampl = noise1;
data12 = ft_freqsimulation(cfg);

data12 = ft_freqsimulation(cfg);



for t = 1:cfg.numtrl,
    data2.trial{t}(idx_same) = data12.trial{t}(idx_same);
end

data1.trialinfo = ones(cfg.numtrl,1);
data2.trialinfo = 2*ones(cfg.numtrl,1);

% cfg = [];
% cfg.trialdef.pre    = 0.5; % s
% cfg.trialdef.post   = 0.5; % s
% cfg.trialdef.numtrl = numtrl;
% cfg.trialdef.type = 1;
% cfg.trialfun = 'test_FT_multiple_comparisons_trialfun';
% 
% [cfgt] = ft_definetrial(cfg);
% data1=ft_redefinetrial(cfgt,data1);
% 
% cfg = [];
% cfg.trialdef.pre    = 0.5; % s
% cfg.trialdef.post   = 0.5; % s
% cfg.trialdef.numtrl = numtrl;
% cfg.trialdef.type = 2;
% cfg.trialfun = 'test_FT_multiple_comparisons_trialfun';
% 
% [cfgt] = ft_definetrial(cfg);
% data2=ft_redefinetrial(cfgt,data2);


figure
subplot(2,1,1);
plot(data1.time{1}, cell2mat(data1.trial'));
title('Raw signal 1');
ylabel('Amplitude');
subplot(2,1,2);
plot(data2.time{1}, cell2mat(data2.trial'));
title('Raw signal 2');
ylabel('Amplitude');
xlabel('Time (s)');

data = ft_appenddata([], data1, data2);

% average across trials, so that mask can be added
cfg           = [];
data1a = ft_timelockanalysis(cfg, data1);
data2a = ft_timelockanalysis(cfg, data2);


cfg           = [];
cfg.method    = 'analytic'; % using a parametric test
cfg.statistic = 'ft_statfun_indepsamplesT'; % using independent samples
cfg.correctm  = 'no'; % no multiple comparisons correction
cfg.alpha     = alpha;
cfg.design    = data.trialinfo; % indicating which trials belong to what category
cfg.ivar      = 1; % indicating that the independent variable is found in first row of cfg.design
stat_t_1 = ft_timelockstatistics(cfg, data);

cfg           = [];
cfg.method           = 'montecarlo'; % use montecarlo to permute the data
cfg.statistic        = 'ft_statfun_indepsamplesT'; % function to use when ...
                                                   % calculating the ...
                                                   % parametric t-values
cfg.alpha     = alpha;
cfg.design    = data.trialinfo; % indicating which trials belong to what category
cfg.ivar      = 1; % indicating that the independent variable is found in first row of cfg.design

cfg.correctm         = 'cluster'; % the correction to use
cfg.clusteralpha     = alpha; % the alpha level used to determine whether or ...
                             % not a channel/time pair can be included in a ...
                             % cluster
cfg.alpha            = alpha/2; % corresponds to an alpha level of 0.05, since ...
                              % two tests are made ...
                              % (negative and positive: 2*0.025=0.05)
cfg.numrandomization = 100;  % number of permutations run
stat_t_2 = ft_timelockstatistics(cfg, data);

cfg               = [];
cfg.baseline = 'no';
cfg.maskparameter = 'mask';
cfg.maskstyle     = 'box';
cfg.channel       = 'mix';
cfg.maskfacealpha = 0.8; % transparency of mask


figure

subplot(3,1,1)
data1a.mask = stat_t_1.mask; % adding stat mask
ft_singleplotER(cfg, data1a, data2a);
hold on
plot([data1a.time(1), data1a.time(end)], [0 0], 'k--') % hor. line
plot(stat_t_1.time,stat_t_1.prob,'k');
plot(stat_t_1.time(stat_t_1.prob<alpha),stat_t_1.prob(stat_t_1.prob<alpha),'r.');
title('no multiple correction');

subplot(3,1,2)
data1a.mask = stat_t_2.mask; % adding stat mask
ft_singleplotER(cfg, data1a, data2a);
hold on
plot([data1a.time(1), data1a.time(end)], [0 0], 'k--') % hor. line
plot(stat_t_2.time,stat_t_2.prob,'k');
plot(stat_t_2.time(stat_t_2.prob<alpha),stat_t_2.prob(stat_t_2.prob<alpha),'r.');
title('cluster');
