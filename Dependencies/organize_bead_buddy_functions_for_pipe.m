%% Generate mapping functions for bead buddy pipeline

% array that will hold our transformation functions
reg_models = cell(max(max(key_tab.SMLM_ch), max(key_tab.bead_ch)));
res_dir = fullfile(bead_dir, 'res');
[fits_found , fits_path] = check_for_fits(res_dir);


% assemble fitfunction object filepath and load into workspace
load(fits_path)



%%
% this is bead ch val for ref scope channel
beadRefChIdx = key_tab.bead_ch( find(key_tab.isReference == 1));

refChannel = key_tab.SMLM_ch( find(key_tab.isReference == 1));

% loop thru all channels
mapError=1;
for i = 1:numel(key_tab(:,"dye"))
    
%     for ith localizable laser line
    if key_tab.isLocalizable(i) == 1 && key_tab.isReference(i) ~= 1

        
        % get the uM FISH idx for the current laser line
        fshIdx = key_tab.SMLM_ch(i);
        

        % get the uM bead idx that is for the same laser line
        beadIdx = key_tab.bead_ch(i);

        % get the reg model that goes from current laser line to the
        % reference laser line (using bead idxs)
        myRegModel = myFits{beadIdx, beadRefChIdx};
        
        % store the appropriate fit function in terms of the FISH indexes
        reg_models{fshIdx, refChannel} = myRegModel;


    end
end
% 
% if mapError == 1
%     error('Mapping failed: Double check your FISH reference channel and channel key')
% end

% Report to user
disp(' ');
disp('Registration functions have been mapped to Loc Data channels');
disp('cell ij contains function mapping Loc Data ch i to ch j (the reference channel')
disp(' ');
disp(reg_models);


%% FUNCTIONS

function [] = miji_process_bead_hyperstacks(bead_path_list, split_ch_dir)

% check if split_ch_dir is already populated
if check_for_file(split_ch_dir, 'tif') 
    disp('You already have stacks in your split channel dir:')
    disp(split_ch_dir)
    disp('Empty the folder in order to proceed')
else

% run miji!
Miji();
MIJ.help %tells u what commands u can do



disp('');
disp('~~~~~~~~~');
disp('Splitting multi-channel images');
disp('');


% loop thru FISH img dirs one condition at a time
for i = 1:numel(bead_path_list)

    my_FOVpath = bead_path_list(i);

    %MIJI open my_FOVpath 

    disp('');

    disp('Processing')
    disp(my_FOVpath);

    % close all open windows
    MIJ.closeAllWindows();

    % open the image in FIJI
    MIJ.run('Open...', strcat('path=', my_FOVpath));

    % get the title
    curTitle = char(MIJ.getCurrentTitle());

    MIJ.run('')

    % split channels
    MIJ.run('Split Channels');

    splitImg_list = MIJ.getListImages;

    nCh = numel(splitImg_list);

    for m=1:nCh


        MIJ.selectWindow(strcat('C', string(m), '-', curTitle));

%             if doFlatField ==1 
%                 MIJ.run('Pseudo flat field correction', 'blurring=50 hide stack');
%             end

        saveTitle = char(MIJ.getCurrentTitle);
        MIJ.run('Save', strcat('Tiff..., path=[', fullfile(split_ch_dir, saveTitle), ']'));
        
    
    end
end
MIJ.closeAllWindows();
MIJ.exit();

end
end


% retruns list of paths to all tifs in input dir
function img_path_list = get_tif_list(my_dir)

all_files = dir(fullfile(my_dir, '*.tif'));

imgs_tab = struct2table(all_files);

img_paths = strings(size(imgs_tab.folder));

for i =1:numel(img_paths)
    
    cur_path = fullfile( string(imgs_tab.folder(i)), string(imgs_tab.name(i)) );
    
    img_paths(i) = string(cur_path);
end
    
img_path_list = img_paths;

end


% returns paths to bead dir and dir of subdir raw images
function [bead_dir, bead_img_dir] = organize_bead_imgs(project_dir)

% lets find all the bead images
all_files = dir(fullfile(project_dir, '**/*.*'));

all_files = struct2table(all_files);

% get sub table with only images
img_idx = contains(all_files.name, 'tif');
bead_idx = contains(all_files.name, 'bead', 'IgnoreCase', 1);
bead_img_table = all_files(img_idx & bead_idx, :);

% for checking redundant file manipulation
d_name = fullfile(project_dir, 'beads');
bead_subdir = d_name;
bead_img_dir = fullfile(d_name, 'raw_imgs');

% check if beads have already been moved
if strcmp(char(bead_img_table.folder(1)), bead_subdir) || strcmp(char(bead_img_table.folder(1)), bead_img_dir)
    disp('~~~~~')
    disp('bead imgs have already been moved to:')
    disp(bead_img_dir )
    
    bead_dir = bead_subdir;
    bead_img_dir = bead_img_dir;
    
    
else

   
    

    % see if beads are already in subfolder of project dir
    if strcmp(char(bead_img_table.folder(1)), project_dir)

        exists_subdir = 0; % there is no subdir for beads 
    else
        exists_subdir = 1; % there is a subdir for beads
        bead_subdir = char(bead_img_table.folder(1)); % save it
        bead_img_dir = fullfile(bead_subdir, 'raw_imgs');
        mkdir(bead_img_dir)
    end

    % if no subdir we will make one
    if exists_subdir == 0 
        d_name = fullfile(project_dir, 'beads');
        mkdir(d_name)
        bead_subdir = d_name;
        bead_img_dir = fullfile(d_name, 'raw_imgs');
        mkdir(bead_img_dir)

    end
    % rename bead files using parent dir
    [~,parent] = fileparts(project_dir); 

    % now systematically rename the FOVS
    % lets find all the bead images
    % get sub table with only images
    for i = 1:numel(bead_img_table.name)

        cur_path = fullfile( string(bead_img_table.folder(i)), string(bead_img_table.name(i)));

        new_path = fullfile( string(bead_img_dir), strcat(parent, '_BEADS', '-', string(i), '.tif') );

        if ~strcmp(cur_path, new_path)

            movefile(cur_path, new_path);
        end

    end

    bead_dir = bead_subdir;
    bead_img_dir = bead_img_dir;

    disp('~~~~~')
    disp('bead image files renamed and organized in directory:')
    disp(bead_subdir)
    disp('')
end
end


function [] = verify_channel_key(key_tab)
% make sure only one ch desiganted as ref
if sum(key_tab.isReference) ~= 1
    disp('Please revise your channel key')
    error('You must designate only one channel as your reference channel')
    
end

% make sure ref channel is localizable
if key_tab.isLocalizable( key_tab.isReference ~= 1 )
    disp('Please revise your channel key')
    error('Your reference channel must be localizable')
end
end

function [fits_found, fits_path] = check_for_fits(project_dir)
% check for fits file
[fits_found, fits_path] = check_for_file(project_dir, 'fits.mat');


% report to user 
if fits_found == 1
    disp('~~~~~')
    disp('fits.mat found at:')
    disp(fits_path)
else
    disp('~~~~~')
    disp('no fits.mat file present in:');
    disp(project_dir)

end
end


function [cfg_found, cfg_path] = check_for_cfg(project_dir)
% check for .ini file
[cfg_found, cfg_path] = check_for_file(project_dir, 'ini');


% report to user 
if cfg_found == 1
    disp('~~~~~')
    disp('cfg file found at:')
    disp(cfg_path)
else
    disp('~~~~~')
    disp('no cfg file present in:');
    disp(project_dir)

end
end


function [key_found, key_path] = check_for_channel_key(project_dir)
% check for channel key file
[key_found, key_path] = check_for_file(project_dir, 'key');


% report to user 
if key_found == 1
    disp('~~~~~')
    disp('channel key found at:')
    disp(key_path)
else
    disp('~~~~~')
    disp('no channel key present in:');
    disp(project_dir)

end
end


function [file_found, file_path] = check_for_file(project_dir, pat)
% CHECK FOR file in input dir (does not open subdirs)
file_found = 0;
file_path = [];
% extract file list
list_dir = dir(project_dir);
% get string array of file names
f_names = string( {list_dir.name} );
folders = string( {list_dir.folder} );

% check if a config file (.ini) is present
for i = 1:numel(f_names)
    if contains(f_names(i), pat, 'IgnoreCase', 1)
        file_found = 1;
        
        file_path = fullfile(folders(i), f_names(i));
        
        
        

    end
    
end

     
end

