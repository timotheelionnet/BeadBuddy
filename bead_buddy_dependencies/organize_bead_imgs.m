
% returns paths to bead dir and dir of subdir raw images
function [bead_dir, bead_img_dir] = organize_bead_imgs(project_dir)

% lets find all the bead images
all_files = dir(fullfile(project_dir, '**/*.*'));

all_files = struct2table(all_files);


% get sub table with only images
img_idx = contains(all_files.name, 'tif');
bead_idx = contains(all_files.name, 'bead', 'IgnoreCase', 1);
bead_img_table = all_files(img_idx & bead_idx, :);

if sum(bead_idx) == 0
    error("You have no image file names with substring 'bead' in them")
end

[d,parent] = (fileparts(bead_img_table.folder(1)));

parent = char(parent);
%%

if strcmp(parent, 'raw_imgs') || strcmp(parent, 'split_channels')

    bead_img_dir = fullfile(d,parent);
    bead_img_dir = char(bead_img_dir);
    bead_dir = char(d);
% check if beads have already been moved
    disp('~~~~~')
    disp('bead imgs have already been moved to:');
    disp(bead_img_dir);

    
    
    
    
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
