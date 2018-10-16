function spike_train = test_FT_simulate_spike_train(firingRate,duration,n_trials)
% see https://praneethnamburi.wordpress.com/2015/02/05/simulating-neural-spike-trains/
% Poisson spike train
% firingRate (Hz)
% duration (s)

if nargin < 3,
	n_trials = 1;
end

dt = 1/1000; % s spike train temporal resolution
nBins = duration*1/dt;
spike_train = rand(n_trials, nBins) < firingRate*dt; 