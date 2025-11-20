

function [] = verify_channel_key(key_tab)
% make sure only one ch desiganted as ref
if sum(key_tab.isReference) ~= 1
    disp('Please revise your channel key')
    error('You must designate only one channel as your reference channel')
    
end

% make sure no localizable channel is set to ch 0
if ismember(0, table2array(key_tab(key_tab.isLocalizable, ["bead_ch","SMLM_ch"])))
    disp('Please revise your channel key')
    error('You must designate a non-zero channel value for each channel to be localized')
    
end

% make sure ref channel is localizable
if key_tab.isLocalizable( find(key_tab.isReference)) ~= 1 
    disp('Please revise your channel key')
    error('Your reference channel must be localizable')
end
end
