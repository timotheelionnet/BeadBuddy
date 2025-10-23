classdef bead_buddy_start_GUI < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                 matlab.ui.Figure
        TopPanel                 matlab.ui.container.Panel
        MiddlePanel              matlab.ui.container.Panel
        BottomPanel              matlab.ui.container.Panel
        ProjectDirTextArea       matlab.ui.control.TextArea
        SelectFolderButton       matlab.ui.control.Button
        NumberInput1             matlab.ui.control.NumericEditField
        NumberInput2             matlab.ui.control.NumericEditField
        OptionCheckBox           matlab.ui.control.CheckBox
        VoxelXYLabel             matlab.ui.control.Label
        VoxelZLabel              matlab.ui.control.Label
        DoneButton               matlab.ui.control.Button
        BeadCorrectCheckBox      matlab.ui.control.CheckBox
        LocateFileButton         matlab.ui.control.Button
        DatasetPathTextArea      matlab.ui.control.TextArea
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: SelectFolderButton
        function SelectFolderButtonPushed(app, event)
            % Prompt user to select a folder
            folder = uigetdir;
            if folder ~= 0
                % Display the selected folder path in the text area
                app.ProjectDirTextArea.Value = folder;
                
                % Save the folder path to workspace as 'project_dir'
                assignin('base', 'project_dir', folder);
            end

            % Bring GUI back on top after file dialog
            if isprop(app.UIFigure, 'AlwaysOnTop')
                drawnow;
                app.UIFigure.AlwaysOnTop = true;
            end
        end

        % Button pushed function: LocateFileButton
        function LocateFileButtonPushed(app, event)
            % Prompt user to select a CSV or XLSX file
            [file, path] = uigetfile({'*.csv;*.xlsx', 'CSV/XLSX Files (*.csv, *.xlsx)'});
            if isequal(file, 0)
                % No file was selected
                app.DatasetPathTextArea.Value = 'No file selected';
            else
                % Display the selected file path
                fullpath = fullfile(path, file);
                app.DatasetPathTextArea.Value = fullpath;
                
                % Save the file path to workspace as 'user_input_data_path'
                assignin('base', 'user_input_data_path', fullpath);
            end

            % Bring GUI back on top after file dialog
            if isprop(app.UIFigure, 'AlwaysOnTop')
                drawnow;
                app.UIFigure.AlwaysOnTop = true;
            end
        end

        % Button pushed function: DoneButton
        function DoneButtonPushed(app, event)
            % Get the values from the input fields
            voxel_xy_value = app.NumberInput1.Value;
            voxel_z_value = app.NumberInput2.Value;
            
            % Save the values to the workspace as 'voxel_xy' and 'voxel_z'
            assignin('base', 'voxel_xy', voxel_xy_value);
            assignin('base', 'voxel_z', voxel_z_value);
            
            % Check the state of the "Run Airlocalize?" checkbox and save runAirlocalize
            if app.OptionCheckBox.Value
                runAirlocalize = 1;
            else
                runAirlocalize = 0;
            end
            % Save runAirlocalize to the workspace
            assignin('base', 'runAirlocalize', runAirlocalize);
            
            % Check the state of the "Bead correct your data?" checkbox and save process_user_data
            if app.BeadCorrectCheckBox.Value
                process_user_data = 1;
            else
                process_user_data = 0;
            end
            % Save process_user_data to the workspace
            assignin('base', 'process_user_data', process_user_data);
            
            % Close the app
            delete(app.UIFigure);
        end
        
        % Value changed function: BeadCorrectCheckBox
        function BeadCorrectCheckBoxValueChanged(app, event)
            % Enable or disable the dataset selection depending on the checkbox state
            if app.BeadCorrectCheckBox.Value
                app.LocateFileButton.Enable = 'on';
                app.DatasetPathTextArea.Enable = 'on';
            else
                app.LocateFileButton.Enable = 'off';
                app.DatasetPathTextArea.Enable = 'off';
            end
        end

    end

    % App initialization and construction
    methods (Access = private)

        % Create and configure components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 400 400];
            app.UIFigure.Name = 'Bead Buddy Start GUI';

            % Create TopPanel for file selection
            app.TopPanel = uipanel(app.UIFigure);
            app.TopPanel.Position = [0 270 400 130];

            % Create ProjectDirTextArea (smaller text area for file path)
            app.ProjectDirTextArea = uitextarea(app.TopPanel);
            app.ProjectDirTextArea.Editable = 'off';
            app.ProjectDirTextArea.Position = [20 70 250 30];
            app.ProjectDirTextArea.Value = 'Select a project folder...';

            % Create SelectFolderButton
            app.SelectFolderButton = uibutton(app.TopPanel, 'push');
            app.SelectFolderButton.ButtonPushedFcn = createCallbackFcn(app, @SelectFolderButtonPushed, true);
            app.SelectFolderButton.Position = [280 70 100 30];
            app.SelectFolderButton.Text = 'Select Folder';

            % Create MiddlePanel for input fields, labels, and checkbox
            app.MiddlePanel = uipanel(app.UIFigure);
            app.MiddlePanel.Position = [0 120 400 150];

            % Create VoxelXYLabel
            app.VoxelXYLabel = uilabel(app.MiddlePanel);
            app.VoxelXYLabel.HorizontalAlignment = 'right';
            app.VoxelXYLabel.Position = [50 90 100 22];
            app.VoxelXYLabel.Text = 'Voxel XY (nm)';

            % Create NumberInput1 (numeric input field for Voxel XY)
            app.NumberInput1 = uieditfield(app.MiddlePanel, 'numeric');
            app.NumberInput1.Position = [160 90 100 22];
            app.NumberInput1.Value = 0;

            % Create VoxelZLabel
            app.VoxelZLabel = uilabel(app.MiddlePanel);
            app.VoxelZLabel.HorizontalAlignment = 'right';
            app.VoxelZLabel.Position = [50 50 100 22];
            app.VoxelZLabel.Text = 'Voxel Z (nm)';

            % Create NumberInput2 (numeric input field for Voxel Z)
            app.NumberInput2 = uieditfield(app.MiddlePanel, 'numeric');
            app.NumberInput2.Position = [160 50 100 22];
            app.NumberInput2.Value = 0;

            % Create OptionCheckBox (checkbox) with default checked state
            app.OptionCheckBox = uicheckbox(app.MiddlePanel);
            app.OptionCheckBox.Text = 'Run Airlocalize?';
            app.OptionCheckBox.Position = [160 20 120 22];
            app.OptionCheckBox.Value = true; % Checked by default

            % Create BottomPanel for "Bead correct your data?" and dataset selection
            app.BottomPanel = uipanel(app.UIFigure);
            app.BottomPanel.Position = [0 0 400 120];

            % Create BeadCorrectCheckBox with default checked state
            app.BeadCorrectCheckBox = uicheckbox(app.BottomPanel);
            app.BeadCorrectCheckBox.Text = 'Bead correct your data?';
            app.BeadCorrectCheckBox.Position = [20 80 150 22];
            app.BeadCorrectCheckBox.Value = true; % Checked by default
            app.BeadCorrectCheckBox.ValueChangedFcn = createCallbackFcn(app, @BeadCorrectCheckBoxValueChanged, true);

            % Create DatasetPathTextArea (text area for dataset file path)
            app.DatasetPathTextArea = uitextarea(app.BottomPanel);
            app.DatasetPathTextArea.Editable = 'off';
            app.DatasetPathTextArea.Position = [20 50 360 30];
            app.DatasetPathTextArea.Value = 'No file selected';

            % Create LocateFileButton for file selection aligned to the left
            app.LocateFileButton = uibutton(app.BottomPanel, 'push');
            app.LocateFileButton.ButtonPushedFcn = createCallbackFcn(app, @LocateFileButtonPushed, true);
            app.LocateFileButton.Position = [20 10 150 30];
            app.LocateFileButton.Text = 'Locate Dataset (csv, xlsx)';

            % Create DoneButton in the bottom right corner
            app.DoneButton = uibutton(app.UIFigure, 'push');
            app.DoneButton.ButtonPushedFcn = createCallbackFcn(app, @DoneButtonPushed, true);
            app.DoneButton.Position = [300 10 70 30];
            app.DoneButton.Text = 'Done';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';

            % In createComponents after making UIFigure
            try
                app.UIFigure.AlwaysOnTop = true; % R2023b+
            catch
                % fallback: do nothing or use JavaFrame trick if available
            end

        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = bead_buddy_start_GUI

            % Create and configure components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Pause script execution until the GUI is closed
            waitfor(app.UIFigure)
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end
