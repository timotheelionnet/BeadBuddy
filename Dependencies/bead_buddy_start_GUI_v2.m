classdef bead_buddy_start_GUI_v2 < matlab.apps.AppBase

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

        % Helper to bring GUI back to front
        function bringToFront(app)
            drawnow;
            try
                if isprop(app.UIFigure, 'AlwaysOnTop')
                    app.UIFigure.AlwaysOnTop = true;
                end
            catch
                % ignore if not supported
            end
            pause(0.05); % let OS process focus change
            try
                figure(app.UIFigure); % request MATLAB focus
            catch
            end
        end

        % Button pushed function: SelectFolderButton
        function SelectFolderButtonPushed(app, event)
            folder = uigetdir;
            if folder ~= 0
                app.ProjectDirTextArea.Value = folder;
                assignin('base', 'project_dir', folder);
            end
            app.bringToFront();
        end

        % Button pushed function: LocateFileButton
        function LocateFileButtonPushed(app, event)
            [file, path] = uigetfile({'*.csv;*.xlsx', 'CSV/XLSX Files (*.csv, *.xlsx)'});
            if isequal(file, 0)
                app.DatasetPathTextArea.Value = 'No file selected';
            else
                fullpath = fullfile(path, file);
                app.DatasetPathTextArea.Value = fullpath;
                assignin('base', 'user_input_data_path', fullpath);
            end
            app.bringToFront();
        end

        % Button pushed function: DoneButton
        function DoneButtonPushed(app, event)
            voxel_xy_value = app.NumberInput1.Value;
            voxel_z_value = app.NumberInput2.Value;
            assignin('base', 'voxel_xy', voxel_xy_value);
            assignin('base', 'voxel_z', voxel_z_value);

            runAirlocalize = double(app.OptionCheckBox.Value);
            assignin('base', 'runAirlocalize', runAirlocalize);

            process_user_data = double(app.BeadCorrectCheckBox.Value);
            assignin('base', 'process_user_data', process_user_data);

            delete(app.UIFigure);
        end

        % Value changed function: BeadCorrectCheckBox
        function BeadCorrectCheckBoxValueChanged(app, event)
            if app.BeadCorrectCheckBox.Value
                app.LocateFileButton.Enable = 'on';
                app.DatasetPathTextArea.Enable = 'on';
            else
                app.LocateFileButton.Enable = 'off';
                app.DatasetPathTextArea.Enable = 'off';
            end
            app.bringToFront();
        end

    end

    % App initialization and construction
    methods (Access = private)

        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 400 400];
            app.UIFigure.Name = 'Bead Buddy Start GUI';

            % Panels and controls
            app.TopPanel = uipanel(app.UIFigure, 'Position', [0 270 400 130]);

            app.ProjectDirTextArea = uitextarea(app.TopPanel, ...
                'Editable', 'off', ...
                'Position', [20 70 250 30], ...
                'Value', 'Select a project folder...');

            app.SelectFolderButton = uibutton(app.TopPanel, 'push', ...
                'ButtonPushedFcn', createCallbackFcn(app, @SelectFolderButtonPushed, true), ...
                'Position', [280 70 100 30], ...
                'Text', 'Select Folder');

            app.MiddlePanel = uipanel(app.UIFigure, 'Position', [0 120 400 150]);

            app.VoxelXYLabel = uilabel(app.MiddlePanel, ...
                'HorizontalAlignment', 'right', ...
                'Position', [50 90 100 22], ...
                'Text', 'Voxel XY (nm)');

            app.NumberInput1 = uieditfield(app.MiddlePanel, 'numeric', ...
                'Position', [160 90 100 22], 'Value', 0);

            app.VoxelZLabel = uilabel(app.MiddlePanel, ...
                'HorizontalAlignment', 'right', ...
                'Position', [50 50 100 22], ...
                'Text', 'Voxel Z (nm)');

            app.NumberInput2 = uieditfield(app.MiddlePanel, 'numeric', ...
                'Position', [160 50 100 22], 'Value', 0);

            app.OptionCheckBox = uicheckbox(app.MiddlePanel, ...
                'Text', 'Run Airlocalize?', ...
                'Position', [160 20 120 22], ...
                'Value', true);

            app.BottomPanel = uipanel(app.UIFigure, 'Position', [0 0 400 120]);

            app.BeadCorrectCheckBox = uicheckbox(app.BottomPanel, ...
                'Text', 'Bead correct your data?', ...
                'Position', [20 80 150 22], ...
                'Value', true, ...
                'ValueChangedFcn', createCallbackFcn(app, @BeadCorrectCheckBoxValueChanged, true));

            app.DatasetPathTextArea = uitextarea(app.BottomPanel, ...
                'Editable', 'off', ...
                'Position', [20 50 360 30], ...
                'Value', 'No file selected');

            app.LocateFileButton = uibutton(app.BottomPanel, 'push', ...
                'ButtonPushedFcn', createCallbackFcn(app, @LocateFileButtonPushed, true), ...
                'Position', [20 10 150 30], ...
                'Text', 'Locate Dataset (csv, xlsx)');

            app.DoneButton = uibutton(app.UIFigure, 'push', ...
                'ButtonPushedFcn', createCallbackFcn(app, @DoneButtonPushed, true), ...
                'Position', [300 10 70 30], ...
                'Text', 'Done');

            % Show the figure
            app.UIFigure.Visible = 'on';
            app.bringToFront();
        end
    end

    % App creation and deletion
    methods (Access = public)

        function app = bead_buddy_start_GUI_v2
            createComponents(app)
            registerApp(app, app.UIFigure)
            waitfor(app.UIFigure)
        end

        function delete(app)
            delete(app.UIFigure)
        end
    end
end
