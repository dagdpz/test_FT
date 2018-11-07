function out = ig_limit_range_min_max(s,mini,maxi)
%IG_LIMIT_RANGE_MIN_MAX		- limit range between min and max
% s can be any matrix

out = arrayfun(@(x) min( maxi, max( mini, x ) ),s);
% out = min( maxi, max( mini, s ) ); % for single element