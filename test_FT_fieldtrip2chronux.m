function croD = test_FT_fieldtrip2chronux(ftD,datatype)
%test_FT_fieldtrip2chronux  - converts fieldtrip data format to chronux
%
% USAGE:
% usage example1;
% usage example2;
%		...
%
% INPUTS:
%		input 1		- exlanation
%		...
%
% OUTPUTS:
%		output1		- explanation
%		...
%
% REQUIRES:	...
%
% See also ...
%
%
% Author(s):	I.Kagan, DAG, DPZ
% URL:		http://www.dpz.eu/dag
%
% Change log:
% yyyymmdd:	Created function (Author(s) firstname familyname)
% ...
% $Revision: 1.0 $  $Date: 2018-10-17 17:03:03 $

% ADDITIONAL INFO:
% ...
%%%%%%%%%%%%%%%%%%%%%%%%%[DAG mfile header version 1]%%%%%%%%%%%%%%%%%%%%%%%%% 

% LFP Fieldtrip format: data.trial
% LFP Chronux format: time x trials
% for single LFP channel
if strcmp(datatype,'lfp'),
	croD = (cell2mat(ftD.trial'))';
	
else % spikes
	spikes = (cell2mat(ftD.trial'))';
	if length(unique(spikes))==2, % spike train: 0s and 1s
		for t = 1:size(ftD.trial,2), % for each trial
			croD(t).times = (ftD.time{t}(ftD.trial{t}==1))'; % spike times should be a column!
		end
	end
end
	