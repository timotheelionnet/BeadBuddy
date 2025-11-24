%% BEADBUDDY

% bead image files must be tifs
% put all bead image tiffs in a subfolder within your project folder. Bead folder name should contain the substring "bead" 
% bead image files must be saved with substring 'bead' in filename
% bead image files cannot have '-' anywere in the name
% user input datatable is expected to have units of pixels
% If 3D, then columns: channel x, y, z, FOV
% If 2D, then columns: channel x, y, FOV


%% Check for required Toolboxes
% List of required toolboxes
requiredToolboxes = { ...
    'Curve Fitting Toolbox',...
    'Image Processing Toolbox',...
    'Optimization Toolbox',...
    'Statistics and Machine Learning Toolbox'

};

% Get installed toolboxes
v = ver;
installedToolboxes = {v.Name};

% Check for missing toolboxes
missingToolboxes = setdiff(requiredToolboxes, installedToolboxes);

if ~isempty(missingToolboxes)
    error(['The following required MATLAB add-ons/toolboxes are missing:' newline ...
        strjoin(missingToolboxes, newline) newline ...
        'Please install them to run this script.']);
end

%% set working directory + add paths
code_dir = fileparts(which('BeadBuddy_MAIN.m'));
cd(code_dir)

% Generate a path string including all subfolders
pathWithSubfolders = genpath(code_dir);

% Add the path to MATLAB's search path
addpath(pathWithSubfolders);

% Display a message indicating that the path has been added
disp(['Added the following directory to the MATLAB path: ' code_dir]);

% See if we have correct FIJI path already saved in our txt file
% FIJI_path_text = fullfile('FIJI_path.txt');
% FIJI_path = fileread(FIJI_path_text);
FIJI_path = fullfile(code_dir, "Dependencies/MIJI/");

FIJI_path = checkAndPromptForFIJIPath(FIJI_path, code_dir);


% use the GUI to get input info
bead_buddy_start_GUI_v3

voxSize = [voxel_xy, voxel_xy, voxel_z];


%~~~RUN~~~%
MIJ_installer


%% CHANNEL KEY
% check for channel_key

[key_found, key_path] = check_for_channel_key(project_dir);


% generate channel key template test GUI
if ~key_found
    % key_path = generate_key_tab(project_dir, 5);
    generate_channel_key_gui_V2(project_dir)
    
    disp('please fill out the table template via GUI')
    [key_found, key_path] = check_for_channel_key(project_dir);
    
end

% load key
key_tab = readtable(key_path);

% check for errors in key
verify_channel_key(key_tab)

% interpret channel key and write a cfg file
[~,parent] = fileparts(project_dir);
ref_ch = key_tab.SMLM_ch( key_tab.isReference == 1);
bead_ref_ch = key_tab.bead_ch( key_tab.isReference == 1);
bead_ch_vals = key_tab.bead_ch( key_tab.bead_ch ~=0);


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
    disp('~~~~~')
    disp("BeadBuddy Complete :)")
    disp('~~~~~')
    return
end

%% Apply to user specified table (maybe have an easy way to apply to image format too)

% % read table with col names: channel, x,y,z or just channel, x,y columns (Assumes
% pixel units of input) 
% Reg is performed in nm
% Saved as pixels

% Use the first img in the split ch dir to get size of the image/camera
% sensor in units of pixels
tmp_img_list = get_tif_list(split_ch_dir);
tmp_img = imread(tmp_img_list(1));
sensor_size = size(tmp_img);


% using GUI
raw_data_path = user_input_data_path;
my_data = readtable(user_input_data_path);
disp(head(my_data));

all_fish_channels = unique(my_data.channel);

% get ref channel name

ref_mask = key_tab.isReference == 1;
ref_ch_name = string(key_tab(ref_mask, :).dye);


% % check headers for 2 or 3D data
col_names = my_data.Properties.VariableNames;

% case insensitive check if there is a z column
is_3D = any(strcmpi('z', col_names));
ncol = numel(my_data.Properties.VariableNames);

if is_3D ==1

    

    if ncol ~= 5 % we need input columns :C, x, y, z, FOV
        error('Input data must have columns in order: channel, x, y, z, FOV')
    end

    % rename columns for downstream use
    my_data.Properties.VariableNames = {'channel', 'x', 'y', 'z', 'FOV'};
end


% output_tab initialize
output_tab = my_data(my_data.channel == ref_ch, :);

% get channels in users dataset
reg_ch_list = setdiff(unique(my_data.channel), ref_ch);

% model output dir
model_out_dir = fullfile(project_dir, 'bead_analysis');
make_new_dir(model_out_dir)

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
        disp('~~~~~')
        disp(strcat('Correcting 3D Dataset channel', string(cur_ch)))
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
        x_max = sensor_size(2) * voxSize(1) ;
        y_min = 1;
        y_max = sensor_size(1) * voxSize(2);
    
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

        % get name of channel from key tab

        ch_mask = key_tab.SMLM_ch == cur_ch;
        cur_reg_ch_name = strcat(string(key_tab.dye(ch_mask)), "-", string(key_tab.names(ch_mask)));

        % now we want a heat map
        mesh_dr_reshape =  reshape(mesh_dr, size(X));

        f = figure;
        imagesc(mesh_dr_reshape)
        cb0 = colorbar;
        cb0_lims = cb0.Limits;
        colormap("Turbo")
        ylabel(cb0, 'absolute displacement (nm)')
        my_title = strcat("BeadBuddy model of absolute chromatic error for Ch", string(cur_ch), "-", cur_reg_ch_name, ' vs Reference Ch-', ref_ch_name);
        title(my_title)
        
        save_dir = fullfile(model_out_dir, 'heatmaps');
        make_new_dir(save_dir)
        save_name = fullfile(save_dir, my_title + '.pdf');

        saveas(f, save_name);
        disp('BeadBuddy chromatic error heat map saved to')
        disp(save_name)
        disp('')

        % we want to also generate pixel scale images for dx dy dz

        % Create a linearly spaced vector for each dimension, evaluate fit func at
        % center of each pixel
        x1 = 0.5*voxSize(1):voxSize(1):sensor_size(2)*voxSize(1) - 0.5*voxSize(1);
        y1 = 0.5*voxSize(2):voxSize(2):sensor_size(1)*voxSize(2) - 0.5*voxSize(2);
    
        % Create the 2D grid of points
        [X1, Y1] = meshgrid(x1, y1);
    
        % evaluate current fit functiosn over mesh
        mesh_dx1 = feval(cur_dx_func, [X1(:) Y1(:)]);
        mesh_dy1 = feval(cur_dy_func, [X1(:) Y1(:)]);
        mesh_dz1 = feval(cur_dz_func, [X1(:) Y1(:)]);

        mesh_dx1_reshape = reshape(mesh_dx1, size(X1));
        mesh_dy1_reshape = reshape(mesh_dy1, size(X1));
        mesh_dz1_reshape = reshape(mesh_dz1, size(X1));
        
        % 3d dataset
        fdx = figure;
        imagesc(mesh_dx1_reshape)
        colormap("winter")
        c = colorbar;
        ylabel(c, 'displacement (nm)')
        title(strcat("dx--",cur_reg_ch_name, ' vs reference ch-', ref_ch_name));
        save_dir = fullfile(model_out_dir, 'displacement_tiffs');
        make_new_dir(save_dir)
        save_name = fullfile(save_dir, strcat("dx--",cur_reg_ch_name) + '.tif');
        % saveas(fdx, save_name);

        save_float_tiff(mesh_dx1_reshape, save_dir, strcat("dx--",cur_reg_ch_name, '.tif'));

        fdy = figure;
        imagesc(mesh_dy1_reshape)
        colormap("spring")
        c = colorbar;
        ylabel(c, 'displacement (nm)')
        title(strcat("dy--",cur_reg_ch_name, ' vs reference ch-', ref_ch_name));

        save_name = fullfile(save_dir, strcat("dy--",cur_reg_ch_name) + '.tif');
        % saveas(fdy, save_name);
        save_float_tiff(mesh_dy1_reshape, save_dir, strcat("dy--",cur_reg_ch_name, '.tif'));

        fdz = figure;
        imagesc(mesh_dz1_reshape)
        colormap("summer")
        c = colorbar;
        ylabel(c, 'displacement (nm)')
        title(strcat("dz--",cur_reg_ch_name, ' vs reference ch-', ref_ch_name));

        save_name = fullfile(save_dir, strcat("dz--",cur_reg_ch_name) + '.tif');
        % saveas(fdz, save_name);

        save_float_tiff(mesh_dz1_reshape, save_dir, strcat("dz--",cur_reg_ch_name) + '.tif');

        disp(' ')
        disp('dx,dy,dz displacement fields (units = nm) as a function of pixel saved as signed tiff images to:')
        disp(save_dir)
        disp('')
        




    else % 2D dataset only correct x and y columns

        if ncol ~= 4 % we need input columns C, x, y, FOV
            error('Input data must have columns in order: channel, x, y, FOV')
        end

        % no z 
        voxSize = voxSize(1:2);

        % apply corrections
        old_x = cur_ch_tab.x;
        old_y = cur_ch_tab.y;
        

        p = [old_x, old_y];

        p_nm = convert_loc_pix_to_nm(p , voxSize);

        dx = feval(cur_dx_func, p_nm(:, 1:2));
        dy = feval(cur_dy_func, p_nm(:, 1:2));
   

        p_nm_new = p_nm - [dx dy];

        p_new = convert_loc_nm_to_pix(p_nm_new, voxSize);

        % lets evaluate the fit function over a 2D grid of test points-----------------------
    
        % Specify the domain
        x_min = 1;
        x_max = sensor_size(2) * voxSize(1) ;
        y_min = 1;
        y_max = sensor_size(1) * voxSize(2);
    
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

        mesh_dr = sqrt(mesh_dx.^2 + mesh_dy.^2);

        % get name of channel from key tab

        ch_mask = key_tab.SMLM_ch == cur_ch;
        cur_reg_ch_name = strcat(string(key_tab.dye(ch_mask)), "-", string(key_tab.names(ch_mask)));
        

        % now we want a heat map
        mesh_dr_reshape =  reshape(mesh_dr, size(X));

        f = figure;
        imagesc(mesh_dr_reshape)
        cb0 = colorbar;
        cb0_lims = cb0.Limits;
        colormap("Turbo")
        ylabel(cb0, 'absolute displacement (nm)')
        my_title = strcat("BeadBuddy model of absolute chromatic error for Ch", string(cur_ch), "-", cur_reg_ch_name, ' vs Reference Ch-', ref_ch_name);
        title(my_title)
        
        save_dir = fullfile(model_out_dir, 'heatmaps');
        make_new_dir(save_dir)
        save_name = fullfile(save_dir, my_title + '.pdf');

        saveas(f, save_name);
        disp('BeadBuddy chromatic error heat map saved to')
        disp(save_name)
        disp('')

        % we want to also generate pixel scale images for dx dy

        % Create a linearly spaced vector for each dimension, evaluate fit func at
        % center of each pixel
        x1 = 0.5*voxSize(1):voxSize(1):sensor_size(2)*voxSize(1) - 0.5*voxSize(1);
        y1 = 0.5*voxSize(2):voxSize(2):sensor_size(1)*voxSize(2) - 0.5*voxSize(2);
    
        % Create the 2D grid of points
        [X1, Y1] = meshgrid(x1, y1);
    
        % evaluate current fit functiosn over mesh
        mesh_dx1 = feval(cur_dx_func, [X1(:) Y1(:)]);
        mesh_dy1 = feval(cur_dy_func, [X1(:) Y1(:)]);
      

        mesh_dx1_reshape = reshape(mesh_dx1, size(X1));
        mesh_dy1_reshape = reshape(mesh_dy1, size(X1));
        
        %2D datset
        fdx = figure;
        imagesc(mesh_dx1_reshape)
        colormap("winter")
        c = colorbar;
        ylabel(c, 'displacement (nm)')
        title(strcat("dx--",cur_reg_ch_name, ' vs reference ch-', ref_ch_name));
        save_dir = fullfile(model_out_dir, 'displacement_tiffs');
        make_new_dir(save_dir)
        save_name = fullfile(save_dir, strcat("dx--",cur_reg_ch_name) + '.tif');
        % saveas(fdx, save_name);

        save_float_tiff(mesh_dx1_reshape, save_dir, strcat("dx--",cur_reg_ch_name, '.tif'));

        fdy = figure;
        imagesc(mesh_dy1_reshape)
        colormap("spring")
        c = colorbar;
        ylabel(c, 'displacement (nm)')
        title(strcat("dy--",cur_reg_ch_name, ' vs reference ch-', ref_ch_name));

        save_name = fullfile(save_dir, strcat("dy--",cur_reg_ch_name) + '.tif');
        % saveas(fdy, save_name);
        save_float_tiff(mesh_dy1_reshape, save_dir, strcat("dy--",cur_reg_ch_name, '.tif'));



        
        disp('dx, dy, tiff images saved to:')
        disp(save_dir)
        disp(' ')
 


    end

 

    
    % generate a subtable for each corrected channel
    cur_ch_output = array2table( [cur_ch_tab.channel, p_new ,cur_ch_tab.FOV]);
    cur_ch_output.Properties.VariableNames = output_tab.Properties.VariableNames;


    % concatenate it to the reference channel data
    output_tab = vertcat(output_tab, cur_ch_output);
    


end
close all


%% PLOTTING AND RESULTS


% color palette for CDFs and other plots. 
color_pal = lines((numel(all_fish_channels) * numel(all_fish_channels))*0.5  + 3);

%% Register user data visualization

[d,f,ext] = fileparts(raw_data_path);

res_dir = fullfile(project_dir, 'results');

make_new_dir(res_dir)

% loop thru raw vs corrected data for quiver plot
for i = 1:numel(reg_ch_list)

    cur_ch = reg_ch_list(i);

    sub_tab_old = my_data(my_data.channel == cur_ch, :);
    
    % handle if user data is 3d or 2d
    if is_3D == 1

        r_old = [sub_tab_old(:,"x"), sub_tab_old(:,"y"), sub_tab_old(:,"z")];
        
    
        sub_tab_new = output_tab(output_tab.channel == cur_ch, :);
    
        r_new = [sub_tab_new(:,"x"), sub_tab_new(:,"y"), sub_tab_new(:,"z")];
        

    else % 2D only

        r_old = [sub_tab_old(:,"x"), sub_tab_old(:,"y")];
        
    
        sub_tab_new = output_tab(output_tab.channel == cur_ch, :);
    
        r_new = [sub_tab_new(:,"x"), sub_tab_new(:,"y")];
        

    end



    % QUIVER xyz analysis for plots in nm

    cur_color = color_pal(i,:);

    xy_old = [r_old.x, r_old.y];

    xy_old_nm = convert_loc_pix_to_nm(xy_old, voxSize(1:2));

    xy_new = [r_new.x, r_new.y];

    xy_new_nm = convert_loc_pix_to_nm(xy_new, voxSize(1:2));

    dxy = xy_new - xy_old;

    dxy_nm = xy_new_nm - xy_old_nm;

    if is_3D == 1

        dz = r_new.z - r_old.z;

        dz_nm = dz * voxSize(3);
    end
    
    % make quiver plot
    f1 = figure;
    % scatter(r_old.x, r_old.y)
    % hold on
    % scatter(r_new.x, r_new.y)

    
    quiver_scale = 2; % svale up quiver size for visualization purposes
    q = quiver(xy_old_nm(:,1), xy_old_nm(:,2), dxy_nm(:,1), dxy_nm(:,2), 2) ;
    q.Color = cur_color;

    % need to set axis limits
    xlim([0, max( xy_old_nm(:,1) + dxy_nm(:,1) )])
    ylim([0, max( xy_old_nm(:,2) + dxy_nm(:,2) )])
    
    myTitle = sprintf('Quiver Ch%d vs Ref Ch%d--XY Correction', [reg_ch_list(i), ref_ch] );
    title(myTitle);
    subtitle( sprintf('Quiver Scale = %d', quiver_scale));

    xlabel('x (nm)')
    ylabel('y (nm)')

    save_dir = fullfile(res_dir, 'registration_viz');
    make_new_dir(save_dir)

    savefig(f1, fullfile(save_dir, myTitle));
    saveas(f1, fullfile(save_dir, strcat(myTitle, '.png') ));

    % dz heat map ------------------------------------
    % if 3d, visualize the the magnitude of z correction
    if is_3D == 1 

        zResids = dz_nm;  % true value is 0, so dz is residuals
    
        %plot residual map
        zResidsMap = min(quantile(zResids,0.9), max(quantile(zResids,0.1),zResids));
        f2 = figure;
        scatter3(xy_new_nm(:,1), xy_new_nm(:,2), zResids, 12*ones(size(r_new.x)),zResidsMap,'filled');
        % quiver(myX, myY, myZ1, myZ2); %plots vectors originating at myX myY and with direction myZ1 myZ2
    
        myTitle = sprintf('Z Correction Ch%d vs Ref Ch%d', [reg_ch_list(i), ref_ch] );
        title(myTitle)
        xlabel('x (nm)')
        ylabel('y (nm)')
        zlabel('z (nm)')
    
    
        savefig(f2, fullfile(save_dir, myTitle));
        saveas(f2, fullfile(save_dir, strcat(myTitle, '.png') ));
    end
  
   
end

disp('~~~~~~~~~~');
disp('Data correction plots saved to:')
disp(save_dir)
disp('~~~~~~~~~~');


close all

%% CDF plots for all combos
% nearest neighbor analysis between channels per FOV, CDF plotted for
% all data

[d,f,ext] = fileparts(raw_data_path);

cdf_dir = fullfile(res_dir, 'NN_analysis_CDF_plts');

if ~exist(cdf_dir, 'dir')
    mkdir(cdf_dir)
end

% in nm
nnThresh = sqrt(2) * max(sensor_size) * voxSize(1);

% only keep mutual nn pairs
mutualOnly = 1;

% CDF props
my_LW = 4;

ctr = 1;
for i = 1:numel(all_fish_channels)
    for j = i+1:numel(all_fish_channels)

        ci = all_fish_channels(i);
        cj = all_fish_channels(j);

        if is_3D == 1
            % get r for ci before and after correction
            sub_tab_old1 = my_data(my_data.channel == ci, :);
            r_old1 = [sub_tab_old1.x, sub_tab_old1.y, sub_tab_old1.z];
            r_old1_nm = convert_loc_pix_to_nm(r_old1, voxSize);
    
    
            sub_tab_new1 = output_tab(output_tab.channel == ci, :);
            r_new1 = [sub_tab_new1.x, sub_tab_new1.y, sub_tab_new1.z];
            r_new1_nm = convert_loc_pix_to_nm(r_new1, voxSize);
    
            % get r for cj before and after correction
            sub_tab_old2 = my_data(my_data.channel == cj, :);
            r_old2 = [sub_tab_old2.x, sub_tab_old2.y, sub_tab_old2.z];
            r_old2_nm = convert_loc_pix_to_nm(r_old2, voxSize);
    
    
            sub_tab_new2 = output_tab(output_tab.channel == cj, :);
            r_new2 = [sub_tab_new2.x, sub_tab_new2.y, sub_tab_new2.z];
            r_new2_nm = convert_loc_pix_to_nm(r_new2, voxSize);
        else % 2D dataset

            % get r for ci before and after correction
            sub_tab_old1 = my_data(my_data.channel == ci, :);
            r_old1 = [sub_tab_old1.x, sub_tab_old1.y];
            r_old1_nm = convert_loc_pix_to_nm(r_old1, voxSize(1:2));
    
    
            sub_tab_new1 = output_tab(output_tab.channel == ci, :);
            r_new1 = [sub_tab_new1.x, sub_tab_new1.y];
            r_new1_nm = convert_loc_pix_to_nm(r_new1, voxSize(1:2));
    
            % get r for cj before and after correction
            sub_tab_old2 = my_data(my_data.channel == cj, :);
            r_old2 = [sub_tab_old2.x, sub_tab_old2.y];
            r_old2_nm = convert_loc_pix_to_nm(r_old2, voxSize(1:2));
    
    
            sub_tab_new2 = output_tab(output_tab.channel == cj, :);
            r_new2 = [sub_tab_new2.x, sub_tab_new2.y];
            r_new2_nm = convert_loc_pix_to_nm(r_new2, voxSize(1:2));
        end


        % compute nn matrix for ci vs cj before and after correction
        
        disp('Before BB Correction')
        nn_res_old = get_nn_matrix_from_2arrays(r_old1_nm, r_old2_nm, nDims, nnThresh, mutualOnly);
        disp('After BB Correction')
        nn_res_new = get_nn_matrix_from_2arrays(r_new1_nm, r_new2_nm, nDims, nnThresh, mutualOnly);
        disp(strcat("NN Threshold =", string(nnThresh), "(nm)"))



        f = figure;

        c1 = cdfplot(nn_res_old(:,1));
        
        hold on
        c2 = cdfplot(nn_res_new(:,1));
        grid off

        c1.LineWidth = my_LW;
        c2.LineWidth = my_LW;

        c1.Color = color_pal(ctr,:) .* [0.4 0.8 0.7];
        c2.Color = color_pal(ctr,:);

        
        myTitle = sprintf('Mutual Nearest Neighbors-User Dataset Ch %d vs Ch %d', [ci, cj]);
        title( myTitle )

        x_text = 'NN Distance (nm)';
        xlabel(x_text)

        legend({'raw', 'corrected'}, Location='east');

        % xlim([0,1])

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

close all
%% save corrected user data

[d,f,ext] = fileparts(raw_data_path);

save_name = fullfile(res_dir, strcat(f,'_BEADCRXN',ext));

writetable(output_tab, save_name);

disp('~~~~~~~~~~');
disp('Corrected user dataset saved in units of pixels to:')
disp(save_name)
disp('~~~~~~~~~~');

     
disp("BeadBuddy Complete :)")
fprintf("See %s for plots and results\n", project_dir);
    


