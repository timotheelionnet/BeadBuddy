function generate_channel_key_gui_V2
    % Create a figure window
    hFig = figure('Position', [100, 100, 700, 400], 'Name', 'Interactive Editable Table GUI', 'CloseRequestFcn', @closeRequestFcn);
    
    % Create a text label for rows
    uicontrol('Style', 'text', ...
              'Position', [50, 330, 120, 30], ...
              'String', 'Total # Ch in Hyperstack', ...
              'FontSize', 10);
          
    % Create an edit box for rows
    hRows = uicontrol('Style', 'edit', ...
                      'Position', [180, 335, 100, 30], ...
                      'FontSize', 10);
                  
    % Create a button to generate table
    hGenerateButton = uicontrol('Style', 'pushbutton', ...
                                'Position', [300, 335, 120, 30], ...
                                'String', 'Generate Table', ...
                                'FontSize', 10, ...
                                'Callback', @generate_table_callback);
          
    % Create a button to save the table
    hSaveButton = uicontrol('Style', 'pushbutton', ...
                            'Position', [300, 300, 120, 30], ...
                            'String', 'Save Table', ...
                            'FontSize', 10, ...
                            'Callback', @save_table_callback, ...
                            'Enable', 'off');
                        
    % Axes for the table
    hAxes = axes('Units', 'pixels', 'Position', [50, 50, 600, 250]);
    axis off;

    % Variable to hold the uitable handle
    hTable = [];
    
    % Pause the execution of the script until the figure is closed
    uiwait(hFig);

    function bring_to_front()
        drawnow;
        figure(hFig);  % bring GUI to front
        try
            % Undocumented Java-based always-on-top trick
            jFig = get(handle(hFig), 'JavaFrame');
            jFig.fFigureClient.getWindow.setAlwaysOnTop(true);
            pause(0.05); % short delay so OS registers it
            jFig.fFigureClient.getWindow.setAlwaysOnTop(false);
        catch
            % Ignore if JavaFrame not available (newer MATLAB versions)
        end
    end    
    
    % Nested function for button callback to generate table
    function generate_table_callback(~, ~)
        % Get the number of rows from user input
        numRows = str2double(hRows.String);
        numCols = 7; % Fixed number of columns
        
        % Validate the input
        if isnan(numRows) || numRows <= 0
            errordlg('Please enter a valid positive integer for the number of rows.', 'Invalid Input');
            return;
        end
        
        % Create data with default values for the table
        data = cell(numRows, numCols);
        for i = 1:numRows
            if i == 1
                data{i, 1} = 'Nucleus';     % names (first row)
            else
                data{i, 1} = 'label';       % names (other rows)
            end
            if i == 1
                data{i, 2} = 0;             % isLocalizable (first row)
            else
                data{i, 2} = 1;             % isLocalizable (other rows)
            end
            data{i, 3} = 0;                 % isReference
            data{i, 4} = i;                 % uM_FISH_ch (starting from 1)
            data{i, 5} = i - 1;             % uM_Bead_ch (starting from 0)
            data{i, 6} = 'scope';           % scope
            % Assigning dye values based on row index
            if i <= 5
                dyeValues = {'dapi', 'cy5', 'a594', 'cy3', 'a488'};
                data{i, 7} = dyeValues{i};
            else
                data{i, 7} = 'my fluor';    % dye (for rows beyond the first five)
            end
            bring_to_front();
        end
        
        % Define specific column names
        colNames = {'names', 'isLocalizable', 'isReference', 'uM_FISH_ch', 'uM_Bead_ch', 'scope', 'dye'};
        
        % Display the table in a uitable
        if isempty(hTable)
            hTable = uitable('Parent', hFig, ...
                             'Data', data, ...
                             'ColumnName', colNames, ...
                             'Position', [50, 50, 600, 250], ...
                             'ColumnEditable', true(1, numCols), ...
                             'ColumnWidth', {70});
        else
            set(hTable, 'Data', data, 'ColumnName', colNames, 'Position', [50, 50, 600, 250]);
        end
        
        % Disable and hide the generate button
        set(hGenerateButton, 'Visible', 'off');
        
        % Enable the save button
        set(hSaveButton, 'Enable', 'on');
    end

    % Nested function for button callback to save table
    function save_table_callback(~, ~)
        % Get the data from the uitable
        tableData = get(hTable, 'Data');
        
        % Create column names
        colNames = get(hTable, 'ColumnName');
        
        % Convert the cell array to a table object
        tableObject = cell2table(tableData, 'VariableNames', colNames);
        
        % Ask user to select a directory for saving the .csv file
        dirName = uigetdir('', 'Select a directory to save the .csv file');
        if dirName == 0
            % User canceled the dialog
            return;
        end
        
        % Generate the .csv file name
        [~, saveDirName] = fileparts(dirName);  % Get the last folder name
        csvFileName = fullfile(dirName, [saveDirName, '_channelKey.csv']);
        
        % Write the table to .csv file
        writetable(tableObject, csvFileName);
        
        % Assign the table object to a variable in the base workspace
        assignin('base', 'key_tab', tableObject);
        
        % Display a message indicating that the table has been saved
        msgbox(['Table data has been saved as ', saveDirName, '_channelKey.csv in the directory: ', dirName], 'Table Saved');

        bring_to_front();
        
        % Close the GUI window
        close(hFig);
        
    end

    % Nested function for close request
    function closeRequestFcn(~, ~)
        % Resume execution of the script
        uiresume(hFig);
        % Delete the figure
        delete(hFig);
    end
end
