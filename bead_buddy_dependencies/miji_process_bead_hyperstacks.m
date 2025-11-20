
function [] = miji_process_bead_hyperstacks(bead_path_list, split_ch_dir)

% check if split_ch_dir is already populated
if check_for_file(split_ch_dir, 'tif') | check_for_file(split_ch_dir, 'tif') 
    disp('~~~~~')
    disp('You already have stacks in your split channel dir:')
    disp(split_ch_dir)
    disp('Delete the split channel folder and reset your file tree to initial structure if you need to repeat MIJI image processing')
    disp('~~~~~')
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

    my_FOVpath = fullfile(bead_path_list(i));



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