%% CODE OVERVIEW
% bead image files must be tifs
% bead image files must be saved with substring 'bead' in name
% cannot have '-' anywere in the name


%% TO DO

% - error catches.
% - cfg file directory level




%% USER INPUT COMMENTED OUT FOR GUI USE

% specify path to folder containing your bead image hyperstacks
% project_dir= '/Users/finneganclark/Desktop/20240910_Bead_Buddy_Final_Code/demo_project';

% give vox size in nm (X,Y,Z)
% voxSize = [100,100,250];

% % user must specify the path to their FIJI install
% FIJI_path = 'C:\Users\Lionnet Lab\Documents\Fiji.app';




% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% end of user input
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




%% set working directory + add paths
code_dir = fileparts(matlab.desktop.editor.getActiveFilename);
cd(code_dir)

% Generate a path string including all subfolders
pathWithSubfolders = genpath(code_dir);

% Add the path to MATLAB's search path
addpath(pathWithSubfolders);

% Display a message indicating that the path has been added
disp(['Added the following directory to the MATLAB path: ' code_dir]);

% See if we have correct FIJI path already saved in our txt file
FIJI_path_text = fullfile('FIJI_path.txt');
FIJI_path = fileread(FIJI_path_text);

FIJI_path = checkAndPromptForFIJIPath(FIJI_path, code_dir);


% use the GUI to get input info
bead_buddy_start_GUI

voxSize = [voxel_xy, voxel_xy, voxel_z];


%~~~RUN~~~%
MIJ_installer


%% CHANNEL KEY
% check for channel_key

[key_found, key_path] = check_for_channel_key(project_dir);


% generate channel key template test GUI
if ~key_found
    % key_path = generate_key_tab(project_dir, 5);
    generate_channel_key_gui_V2
    
    disp('please fill out the table template via GUI')
    [key_found, key_path] = check_for_channel_key(project_dir);
    
end

% load key
key_tab = readtable(key_path);

% check for errors in key
verify_channel_key(key_tab)

% interpret channel key and write a cfg file
[~,parent] = fileparts(project_dir);
ref_ch = key_tab.uM_FISH_ch( key_tab.isReference == 1);
bead_ref_ch = key_tab.uM_Bead_ch( key_tab.isReference == 1);
bead_ch_vals = key_tab.uM_Bead_ch( key_tab.uM_Bead_ch ~=0);


key_tab 



%% organize and rename bead images
[bead_dir, bead_img_dir] = organize_bead_imgs(project_dir);


%% WRITE CFG

bead_paths = get_tif_list(bead_img_dir);

% make split channel dir
split_ch_dir = fullfile(bead_dir, 'split_channels');
if ~exist (split_ch_dir, 'dir')
    mkdir(split_ch_dir)
    
    disp('~~~~~')
    disp('bead imgs will be split and saved to:')
    disp(split_ch_dir)
end

% check for config file
% [cfg_found, cfg_path]  = check_for_cfg(project_dir);

% write a cfg file 
cfg_path = bead_cfg_writer_for_pipe(code_dir, split_ch_dir, parent, bead_ch_vals, bead_ref_ch);


%% MIJI

%~~~RUN~~~%
miji_process_bead_hyperstacks(bead_paths, split_ch_dir);


%% AIRLOCALIZE

if runAirlocalize == 1

    %~~~RUN~~~%
    MODULAR_Bead_Analysis_Pipeline_V6 
end

%% Organize bead correction functions


if process_user_data == 1
    %~~~RUN~~~%
    organize_bead_buddy_functions_for_pipe

else
    return
end

%% Apply to user specified table (maybe have an easy way to apply to image format too)

% % read table with col names: channel, x,y,z or just channel, x,y columns (Assumes
% pixel units of input) 
% Reg is performed in nm
% Saved as pixels

% raw_data_path = '/Users/finneganclark/Desktop/20240910_Bead_Buddy_Final_Code/demo_project/bead_buddy_test_data_4ch.xlsx';
% my_data = readtable(raw_data_path);

% using GUI
raw_data_path = user_input_data_path;
my_data = readtable(user_input_data_path);
disp(head(my_data));

all_fish_channels = unique(my_data.channel);


% % check headers for 2 or 3D data
col_names = my_data.Properties.VariableNames;

% case insensitive check if there is a z column
is_3D = any(strcmpi('z', col_names));


if is_3D ==1

    ncol = numel(my_data.Properties.VariableNames);

    if ncol ~= 5 % we need input columns :C, x, y, z, FOV
        error('Input data must have columns in order: channel, x, y, z, FOV')
    end

    % rename columns for downstream use
    my_data.Properties.VariableNames = {'channel', 'x', 'y', 'z', 'FOV'}
end


% output_tab initialize
output_tab = my_data(my_data.channel == ref_ch, :);

% get channels in users dataset
reg_ch_list = setdiff(unique(my_data.channel), ref_ch);
% loop thru
for i = 1:numel(reg_ch_list)

    % extract appropriate reg models 
    cur_ch = reg_ch_list(i);
    cur_reg_model = reg_models(cur_ch, ref_ch);
    cur_reg_model = cur_reg_model{1};
    
    % get the fit objects
    cur_dx_func = cur_reg_model{1}{1};
    cur_dy_func = cur_reg_model{2}{1};

    if is_3D
        cur_dz_func = cur_reg_model{3}{1};

        nDims = 3;
    else
        nDims = 2;
    end


    % subset table for current channel
    cur_ch_tab = my_data( my_data.channel == cur_ch, :);

   
    if is_3D
        % apply corrections
        old_x = cur_ch_tab.x;
        old_y = cur_ch_tab.y;
        old_z = cur_ch_tab.z;

        p = [old_x, old_y, old_z];

        p_nm = convert_loc_pix_to_nm(p , voxSize);

        dx = feval(cur_dx_func, p_nm(:, 1:2));
        dy = feval(cur_dy_func, p_nm(:, 1:2));
        dz = feval(cur_dz_func, p_nm(:, 1:2));

        p_nm_new = p_nm - [dx dy dz];

        p_new = convert_loc_nm_to_pix(p_nm_new, voxSize);


        % lets evaluate the fit function over a grid of test points-----------------------
    
        % Specify the domain
        x_min = 1;
        x_max = 2048 * voxSize(1) ;
        y_min = 1;
        y_max = 2048 * voxSize(2);
    
        % Specify the resolution of the grid
        ds_factor = 128 * voxSize(1); 
        num_points = floor(max(x_max, y_max) / ds_factor);  % Number of points along each dimension
    
        % Create a linearly spaced vector for each dimension
        x = linspace(x_min, x_max, num_points);
        y = linspace(y_min, y_max, num_points);
    
        % Create the 2D grid of points
        [X, Y] = meshgrid(x, y);
    
        % evaluate current fit functiosn over mesh
        mesh_dx = feval(cur_dx_func, [X(:) Y(:)]);
        mesh_dy = feval(cur_dy_func, [X(:) Y(:)]);
        mesh_dz = feval(cur_dz_func, [X(:) Y(:)]);

        mesh_dr = sqrt(mesh_dx.^2 + mesh_dy.^2 + mesh_dz.^2 );

        % now we want a heat map
        mesh_dr_reshape =  reshape(mesh_dr, size(X));

        f = figure;
        imagesc(mesh_dr_reshape)
        cb0 = colorbar;
        cb0_lims = cb0.Limits;
        colormap("Turbo")
        my_title = "Predicted residuals for reg channel " + string(cur_ch);
        title(my_title)
        
        save_dir = fullfile(project_dir, 'Model_visualization');
        make_new_dir(save_dir)
        save_name = fullfile(save_dir, my_title + '.pdf');

        saveas(f, save_name);
        disp('Predicted residuals heat map saved to')
        disp(save_name)

    else % only correct x and y columns

        if ncol ~= 4 % we need input columns :C, x, y, FOV
            error('Input data must have columns in order: channel, x, y, FOV')
        end


        voxSize =voxSize(1:2);

        % apply corrections
        old_x = cur_ch_tab.x;
        old_y = cur_ch_tab.y;
        

        p = [old_x, old_y];

        p_nm = convert_loc_pix_to_nm(p , voxSize);

        dx = feval(cur_dx_func, p_nm(:, 1:2));
        dy = feval(cur_dy_func, p_nm(:, 1:2));
   

        p_nm_new = p_nm - [dx dy];

        p_new = convert_loc_nm_to_pix(p_nm_new, voxSize);
 


    end

 

    
    % generate a subtable for each corrected channel
    cur_ch_output = array2table( [cur_ch_tab.channel, p_new ,cur_ch_tab.FOV]);
    cur_ch_output.Properties.VariableNames = output_tab.Properties.VariableNames;


    % concatenate it to the reference channel data
    output_tab = vertcat(output_tab, cur_ch_output);
    


end

%% PLOTTING AND RESULTS


% color palette for CDFs and other plots. 
color_pal = lines((numel(all_fish_channels) * numel(all_fish_channels))*0.5  + 3);

%% Registration of user data visualization

[d,f,ext] = fileparts(raw_data_path);

residuals_dir = fullfile(d, 'residaul_plots_res');

if ~exist(residuals_dir, 'dir')
    mkdir(residuals_dir)
end

% loop thru raw vs corrected data for quiver plot
for i = 1:numel(reg_ch_list)

    cur_ch = reg_ch_list(i);

    sub_tab_old = my_data(my_data.channel == cur_ch, :);
    
    % handle if user data is 3d or 2d
    if is_3D == 1

        r_old = [sub_tab_old(:,"x"), sub_tab_old(:,"y"), sub_tab_old(:,"z")];
    
        sub_tab_new = output_tab(output_tab.channel == cur_ch, :);
    
        r_new = [sub_tab_new(:,"x"), sub_tab_new(:,"y"), sub_tab_new(:,"z")];

    else

        r_old = [sub_tab_old(:,"x"), sub_tab_old(:,"y")];
    
        sub_tab_new = output_tab(output_tab.channel == cur_ch, :);
    
        r_new = [sub_tab_new(:,"x"), sub_tab_new(:,"y")];

    end



    % xyz analysis for plots

    cur_color = color_pal(i,:);

    xy_old = [r_old.x, r_old.y];

    xy_new = [r_new.x, r_new.y];

    dxy = xy_new - xy_old;

    if is_3D == 1

        dz = r_new.z - r_old.z;
    end
    
    % make quiver plot
    f1 = figure;
    % scatter(r_old.x, r_old.y)
    % hold on
    % scatter(r_new.x, r_new.y)

    
    quiver_scale = 2; % svale up quiver size for visualization purposes
    q = quiver(xy_old(:,1), xy_old(:,2), dxy(:,1), dxy(:,2), 2) ;
    q.Color = cur_color;

    % need to set axis limits
    xlim([0, max( xy_old(:,1) + dxy(:,1) )])
    ylim([0, max( xy_old(:,2) + dxy(:,2) )])
    
    myTitle = sprintf('XY Correction of Channel %d Using Reference Channel %d', [reg_ch_list(i), ref_ch] );
    title(myTitle);
    subtitle( sprintf('Quiver Scale = %d', quiver_scale));

    xlabel('x (pix)')
    ylabel('y (pix)')

    savefig(f1, fullfile(residuals_dir, myTitle));
    saveas(f1, fullfile(residuals_dir, strcat(myTitle, '.png') ));

    % dz heat map ------------------------------------
    % if 3d, visualize the the magnitude of z correction
    if is_3D == 1 

    zResids = dz;  % true value is 0, so dz is residuals
    
        %plot residual map
        zResidsMap = min(quantile(zResids,0.9), max(quantile(zResids,0.1),zResids));
        f2 = figure;
        scatter3(r_new.x, r_new.y, zResids, 12*ones(size(r_new.x)),zResidsMap,'filled');
        % quiver(myX, myY, myZ1, myZ2); %plots vectors originating at myX myY and with direction myZ1 myZ2
    
        myTitle = sprintf('Z Correction of Channel %d Using Reference Channel %d', [reg_ch_list(i), ref_ch] );
        title(myTitle)
        xlabel('x (pix)')
        ylabel('y (pix)')
        zlabel('z (pix')
    
    
        savefig(f2, fullfile(residuals_dir, myTitle));
        saveas(f2, fullfile(residuals_dir, strcat(myTitle, '.png') ));
    end
  
   
end

disp('~~~~~~~~~~');
disp('Registration analsysis plots saved to:')
disp(residuals_dir)
disp('~~~~~~~~~~');




%% CDF plots for all combos
% nearest neighbor analysis between channels per FOV, CDF plotted for
% all data

[d,f,ext] = fileparts(raw_data_path);

cdf_dir = fullfile(d, 'CDF_res');

if ~exist(cdf_dir, 'dir')
    mkdir(cdf_dir)
end

% in pixels
nnThresh = 10;

% only keep mutual nn pairs
mutualOnly = 1;

% CDF props
my_LW = 4;

ctr = 1;
for i = 1:numel(all_fish_channels)
    for j = i+1:numel(all_fish_channels)

        ci = all_fish_channels(i);
        cj = all_fish_channels(j);

        
        % get r for ci before and after correction
        sub_tab_old1 = my_data(my_data.channel == ci, :);
        r_old1 = [sub_tab_old1.x, sub_tab_old1.y, sub_tab_old1.z];


        sub_tab_new1 = output_tab(output_tab.channel == ci, :);
        r_new1 = [sub_tab_new1.x, sub_tab_new1.y, sub_tab_new1.z];

        % get r for cj before and after correction
        sub_tab_old2 = my_data(my_data.channel == cj, :);
        r_old2 = [sub_tab_old2.x, sub_tab_old2.y, sub_tab_old2.z];


        sub_tab_new2 = output_tab(output_tab.channel == cj, :);
        r_new2 = [sub_tab_new2.x, sub_tab_new2.y, sub_tab_new2.z];

        % compute nn matrix for ci vs cj before and after correction

        nn_res_old = get_nn_matrix_from_2arrays(r_old1, r_old2, nDims, nnThresh, mutualOnly);

        nn_res_new = get_nn_matrix_from_2arrays(r_new1, r_new2, nDims, nnThresh, mutualOnly);



        f = figure;

        c1 = cdfplot(nn_res_old(:,1));
        
        hold on
        c2 = cdfplot(nn_res_new(:,1));
        grid off

        c1.LineWidth = my_LW;
        c2.LineWidth = my_LW;

        c1.Color = color_pal(ctr,:) .* [0.4 0.8 0.7];
        c2.Color = color_pal(ctr,:);

        
        myTitle = sprintf('3D Nearest Neighbors Ch %d vs Ch %d', [ci, cj]);
        title( myTitle )

        x_text = '3D Distance (pixels)';
        xlabel(x_text)

        legend({'raw', 'corrected'}, Location='east');

        savefig(f, fullfile(cdf_dir, myTitle));
        saveas(f, fullfile(cdf_dir, strcat(myTitle, '.png') ));


        cdf_table1 = array2table([c1.XData', c1.YData'], "VariableNames",["X", "F(X)"]);
        cdf_table2 = array2table([c2.XData', c2.YData'], "VariableNames",["X", "F(X)"]);

        raw_cdf_fname = fullfile(cdf_dir, sprintf('raw_nn_cdf_ch%d_vs_ch%d.csv', [ci, cj] ) );
        crxn_cdf_fname = fullfile(cdf_dir, sprintf('crxn_nn_cdf_ch%d_vs_ch%d.csv', [ci, cj] ) );

        writetable(cdf_table1, raw_cdf_fname);
        writetable(cdf_table2, crxn_cdf_fname);
        
        
        % color ctr
        ctr = ctr +1 ;
        
    end
end

disp('~~~~~~~~~~');
disp('CDFs saved to:')
disp(cdf_dir)
disp('~~~~~~~~~~');


%% save

[d,f,ext] = fileparts(raw_data_path);

save_name = fullfile(d, strcat(f,'_BEADCRXN',ext));

writetable(output_tab, save_name);

disp('~~~~~~~~~~');
disp('Corrected data saved to:')
disp(save_name)
disp('~~~~~~~~~~');


%% FUNCTIONS
% 
% function [] = miji_process_bead_hyperstacks(bead_path_list, split_ch_dir)
% 
% % check if split_ch_dir is already populated
% if check_for_file(split_ch_dir, 'tif') 
%     disp('~~~~~')
%     disp('You already have stacks in your split channel dir:')
%     disp(split_ch_dir)
%     disp('Delete the split channel folder if you need to MIJI image processing')
%     disp('~~~~~')
% else
% 
% % run miji!
% Miji();
% MIJ.help %tells u what commands u can do
% 
% 
% 
% disp('');
% disp('~~~~~~~~~');
% disp('Splitting multi-channel images');
% disp('');
% 
% 
% % loop thru FISH img dirs one condition at a time
% for i = 1:numel(bead_path_list)
% 
%     my_FOVpath = bead_path_list(i);
% 
%     %MIJI open my_FOVpath 
% 
%     disp('');
% 
%     disp('Processing')
%     disp(my_FOVpath);
% 
%     % close all open windows
%     MIJ.closeAllWindows();
% 
%     % open the image in FIJI
%     MIJ.run('Open...', strcat('path=', my_FOVpath));
% 
%     % get the title
%     curTitle = char(MIJ.getCurrentTitle());
% 
%     MIJ.run('')
% 
%     % split channels
%     MIJ.run('Split Channels');
% 
%     splitImg_list = MIJ.getListImages;
% 
%     nCh = numel(splitImg_list);
% 
%     for m=1:nCh
% 
% 
%         MIJ.selectWindow(strcat('C', string(m), '-', curTitle));
% 
% %             if doFlatField ==1 
% %                 MIJ.run('Pseudo flat field correction', 'blurring=50 hide stack');
% %             end
% 
%         saveTitle = char(MIJ.getCurrentTitle);
%         MIJ.run('Save', strcat('Tiff..., path=[', fullfile(split_ch_dir, saveTitle), ']'));
% 
% 
%     end
% end
% MIJ.closeAllWindows();
% MIJ.exit();
% 
% end
% end
% 
% 
% % retruns list of paths to all tifs in input dir
% function img_path_list = get_tif_list(my_dir)
% 
% all_files = dir(fullfile(my_dir, '*.tif'));
% 
% imgs_tab = struct2table(all_files);
% 
% img_paths = strings(size(imgs_tab.folder));
% 
% for i =1:numel(img_paths)
% 
%     cur_path = fullfile( string(imgs_tab.folder(i)), string(imgs_tab.name(i)) );
% 
%     img_paths(i) = string(cur_path);
% end
% 
% img_path_list = img_paths;
% 
% end
% 
% 
% % returns paths to bead dir and dir of subdir raw images
% function [bead_dir, bead_img_dir] = organize_bead_imgs(project_dir)
% 
% % lets find all the bead images
% all_files = dir(fullfile(project_dir, '**/*.*'));
% 
% all_files = struct2table(all_files);
% 
% % get sub table with only images
% img_idx = contains(all_files.name, 'tif');
% bead_idx = contains(all_files.name, 'bead', 'IgnoreCase', 1);
% bead_img_table = all_files(img_idx & bead_idx, :);
% 
% [d,parent] = (fileparts(bead_img_table.folder(1)));
% 
% parent = char(parent);
% 
% 
% if strcmp(parent, 'raw_imgs') || strcmp(parent, 'split_channels')
% 
%     bead_img_dir = fullfile(d,parent);
%     bead_img_dir = char(bead_img_dir);
%     bead_dir = char(d);
% % check if beads have already been moved
%     disp('~~~~~')
%     disp('bead imgs have already been moved to:');
%     disp(bead_img_dir);
% 
% 
% 
% 
% 
% else
% 
% 
% 
% 
%     % see if beads are already in subfolder of project dir
%     if strcmp(char(bead_img_table.folder(1)), project_dir)
% 
%         exists_subdir = 0; % there is no subdir for beads 
%     else
%         exists_subdir = 1; % there is a subdir for beads
%         bead_subdir = char(bead_img_table.folder(1)); % save it
%         bead_img_dir = fullfile(bead_subdir, 'raw_imgs');
%         mkdir(bead_img_dir)
%     end
% 
%     % if no subdir we will make one
%     if exists_subdir == 0 
%         d_name = fullfile(project_dir, 'beads');
%         mkdir(d_name)
%         bead_subdir = d_name;
%         bead_img_dir = fullfile(d_name, 'raw_imgs');
%         mkdir(bead_img_dir)
% 
%     end
%     % rename bead files using parent dir
%     [~,parent] = fileparts(project_dir); 
% 
%     % now systematically rename the FOVS
%     % lets find all the bead images
%     % get sub table with only images
%     for i = 1:numel(bead_img_table.name)
% 
%         cur_path = fullfile( string(bead_img_table.folder(i)), string(bead_img_table.name(i)));
% 
%         new_path = fullfile( string(bead_img_dir), strcat(parent, '_BEADS', '-', string(i), '.tif') );
% 
%         if ~strcmp(cur_path, new_path)
% 
%             movefile(cur_path, new_path);
%         end
% 
%     end
% 
%     bead_dir = bead_subdir;
%     bead_img_dir = bead_img_dir;
% 
%     disp('~~~~~')
%     disp('bead image files renamed and organized in directory:')
%     disp(bead_subdir)
%     disp('')
% end
% end
% 
% 
% function [] = verify_channel_key(key_tab)
% % make sure only one ch desiganted as ref
% if sum(key_tab.isReference) ~= 1
%     disp('Please revise your channel key')
%     error('You must designate only one channel as your reference channel')
% 
% end
% 
% % make sure ref channel is localizable
% if key_tab.isLocalizable( find(key_tab.isReference)) ~= 1 
%     disp('Please revise your channel key')
%     error('Your reference channel must be localizable')
% end
% end
% 
% function [fits_found, fits_path] = check_for_fits(project_dir)
% % check for fits file
% [fits_found, fits_path] = check_for_file(project_dir, 'fits.mat');
% 
% 
% % report to user 
% if fits_found == 1
%     disp('~~~~~')
%     disp('fits.mat found at:')
%     disp(fits_path)
% else
%     disp('~~~~~')
%     disp('no fits.mat file present in:');
%     disp(project_dir)
% 
% end
% end
% 
% 
% function [cfg_found, cfg_path] = check_for_cfg(project_dir)
% % check for .ini file
% [cfg_found, cfg_path] = check_for_file(project_dir, 'ini');
% 
% 
% % report to user 
% if cfg_found == 1
%     disp('~~~~~')
%     disp('cfg file found at:')
%     disp(cfg_path)
% else
%     disp('~~~~~')
%     disp('no cfg file present in:');
%     disp(project_dir)
% 
% end
% end
% 
% 
% function [key_found, key_path] = check_for_channel_key(project_dir)
% % check for channel key file
% [key_found, key_path] = check_for_file(project_dir, 'key');
% 
% 
% % report to user 
% if key_found == 1
%     disp('~~~~~')
%     disp('channel key found at:')
%     disp(key_path)
% else
%     disp('~~~~~')
%     disp('no channel key present in:');
%     disp(project_dir)
% 
% end
% end
% 
% 
% function [file_found, file_path] = check_for_file(project_dir, pat)
% % CHECK FOR file in input dir (does not open subdirs)
% file_found = 0;
% file_path = [];
% % extract file list
% list_dir = dir(project_dir);
% % get string array of file names
% f_names = string( {list_dir.name} );
% folders = string( {list_dir.folder} );
% 
% % check if a config file (.ini) is present
% for i = 1:numel(f_names)
%     if contains(f_names(i), pat, 'IgnoreCase', 1)
%         file_found = 1;
% 
%         file_path = fullfile(folders(i), f_names(i));
% 
% 
% 
% 
%     end
% 
% end
% 
% 
% end
% 
