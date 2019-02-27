addpath('Y:\Sources\fieldtrip-20190130'); % or other path to fieldtrip
ft_defaults

% http://www.fieldtriptoolbox.org/faq/can_i_prevent_external_toolboxes_from_being_added_to_my_matlab_path/
% (fixed version from above)
[ftver, ftpath] = ft_version;
rmpath(fullfile(ftpath, 'external', 'signal'))
rmpath(fullfile(ftpath, 'external', 'stats'))
rmpath(fullfile(ftpath, 'external', 'images'))

% addpath(genpath('Y:\Sources\chronux_2_12')) % chronux with its subfolders