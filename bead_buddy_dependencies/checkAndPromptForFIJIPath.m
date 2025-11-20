%% Written with chat GPT 
% Finn 6/24/24

function dirPath = checkAndPromptForFIJIPath(dirName, code_dir)

    % Check if the specified directory exists
    if exist(dirName, 'dir') == 7
        fprintf('Directory "%s" found.\n', dirName);
        dirPath = dirName;
    else
        disp('Locate the directory where FIJI is installed');
        
        % Create the GUI to prompt the user to locate the directory
        dirPath = uigetdir(pwd, 'Locate the directory where FIJI is installed');

        % dirPath = fullfile(dirPath, 'Fiji');
        
        % Check if the user cancelled the GUI
        if isequal(dirPath, 0)
            error('Directory not located. Operation cancelled.');
        else
            fprintf('Directory located: %s\n', dirPath);
        end
    end

    % Save the directory path to a text file in the current directory

    filePath = fullfile(code_dir, 'FIJI_path.txt')
    fileID = fopen(filePath, 'wt');  % Open file for writing ('wt' clears existing content)
    if fileID == -1
        error('Could not open FIJI_path.txt for writing.');
    end
    fprintf(fileID, '%s', dirPath);
    fclose(fileID);

    % Save the directory path to the workspace
    assignin('base', 'FIJI_path', dirPath);
end
