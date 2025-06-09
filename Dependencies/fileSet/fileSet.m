classdef fileSet < handle
% v1: changed the formating of the ini file for simplicity

    properties (GetAccess = 'public', SetAccess = 'public')
        
        % file paths
        iniFileName; % the path to the .ini format config file
        outFolder; % the path to the folder where the fileSet is saved. 
        
        % list of file types and conditions in the dataset
        fileTypes; % e.g. {'Loc', 'Mask', 'RawImg'}
        conditions; % e.g. {'FOV', 't', 'Channel'}
        
        % patterns of each file type name
        % 1 x <# of fileTypes> table; each variable has the name of the
        % corresponding fileType; each entry should be a string with the 
        % filename pattern where conditions are encoded within curly
        % brackets:
        % 'Loc' | 'Mask' | 'RawImg'
        % 'xxx/{Channel}-img_t{t}_f{FOV}.loc3' | 'xxx/C1-img_t{t}_f{FOV}_masks.tif' | 'xxx/{Channel}-img_t{t}_f{FOV}.tif'
        patternSet;        
        
        % loading and processing options
        % structure with fields: 
        % recursive (0/1), verbose (0/1), useTerminal (0/1), batchMode
        % (0/1)
        options;
        
        %list of files in the dataset
        % table with following columns:
        % condition_1 | ... | condition_n | fileType_1 ... | fileType_m

        fList;
    end
    
    properties (GetAccess = 'private', SetAccess = 'private')
        defaultConditionNameSingleImg = 'singleFile';
        defaultConditionValueSingleImg = 1;
        defaultOutFileName = 'fileSetList.txt';
        defaultDelimiter = '\t';
        
        % default values for processing options - set in initDefaultOptions
        defaultOptions;
        
        voidEntryString = '<missing>';
        voidEntryValue = NaN;
    end
    
    methods
        %% initialize object
        function obj = fileSet()
            obj.initDefaultOptions;
            obj.initOptions;
        end
        
        %% construct file set from Config file
        function buildFileSetFromConfig(obj,inputFileName)
            if ~isempty(inputFileName)
                obj.iniFileName = inputFileName;
            end
            disp(['Building fileSet from file ',obj.iniFileName,' ...']);
            
            % load output folder info from Config file
            obj.loadOutFolderInfoFromConfig;
            
            % load options from Config file
            obj.loadOptionsFromConfig;
            
            % parse list of conditions from ini config file
            obj.loadConditionsFromConfig;
            
            % load file processing options from config file
            obj.loadFileTypesFromConfig;
            
            % load pattern of paths/names for the different file types and condition from config file
            obj.loadPatternSetFromConfig;
            
            % use input pattern to compile a file list 
            obj.assembleFileListFromPattern;
            
            disp(['Done building fileSet from file ',obj.iniFileName,' ...']);
            disp(' ');
        end
        
        %% create a fileSet object from a table f
        % outFolder can be left empty if not needed.
        % condition argument should contain cell array with names of 
        % the table columns that should be assigned as conditions 
        % properties of the object, e.g. {'FOV','Treatment'}
        % the rest of the columns of the table are assigned as fileTypes.
        % if conditions argument is left empty, a command line prompt will
        % request them from the user.
        function buildFileSetFromFList(obj,f,conditions,outFolder)
            
            obj.fList = f;
            if isempty(conditions)
                disp('Input table columns need to be separated into Conditions and FileTypes.');
                disp('Column names in the input table are: ');
                for i=1:numel(f.Properties.VariableNames)
                    disp([num2str(i),': ',f.Properties.VariableNames{i}]);
                end
                idx = input('Enter which column numbers are Conditions (using array format, e.g. [1,2]): ');
                conditions = f.Properties.VariableNames(idx);
            else
                % make sure all requested conditions exist
                idxToRemove = false(size(conditions));
                for i=1:numel(conditions)
                    if ~ismember(conditions{i},f.Properties.VariableNames)
                        disp(['Warning: there is no column named ',...
                            conditions{i},' in the input table. Cannot use this condition.']);
                        idxToRemove(i) = 1;
                    end
                end
                conditions = conditions(~idxToRemove);
            end
            
            obj.conditions = reshape(conditions,1,numel(conditions));
            obj.fileTypes = f.Properties.VariableNames(...
                ~ismember(f.Properties.VariableNames,conditions));
            obj.fileTypes = reshape(obj.fileTypes,1,numel(obj.fileTypes));
            
            if ~isempty(outFolder)
                obj.outFolder = outFolder;
            end
        end
        
        %% create a fileSet object from a table
        % outFolder can be left empty if not needed.
        % if conditions is left empty, they need to be entered manually 
        % upon command line prompt.
        function buildFileSetFromTableInFile(obj,fileName,conditions,outFolder)
            f = readtable(fileName,'FileType','text',...
                'Delimiter','\t',...
                'ReadVariableNames',true);
            obj.buildFileSetFromFList(f,conditions,outFolder);
        end
        
        %% save file list as a table in the default location 
        % or a specific location (outFname, full path).
        % optional argument: delimiter. (default is tab-delimited; enter ',')
        function saveFileList(obj,outFname,varargin)
            
            if numel(varargin) > 0
                delimiter = varargin{1};
            else
                delimiter = obj.defaultDelimiter;
            end
            if isempty (outFname)
                if isempty(obj.outFolder)
                    disp(['Warning! No output folder specified, ',...
                        'file list will be saved in current working directory: ',...
                        pwd]);
                    disp(['You can specify output folder by setting ',...
                        'the property obj.outFolder = ''path/to/your/folder''']);
                end
                if ~exist(obj.outFolder,'dir')
                    mkdir(obj.outFolder);
                end
                outFname = fullfile(obj.outFolder,obj.defaultOutFileName);
            else
                f = fileparts(outFname);
                if ~isempty(f)
                    if ~exist(f,'dir')
                            mkdir(f);
                    end
                end  
            end
            
            writetable(obj.fList,outFname,'Delimiter',delimiter);
            disp(['Saved file list as text file in: ',outFname]);
        end
        
        %% get the path to the output folder
        function of = getOutFolder(obj)
            of = obj.outFolder;
        end
        
        %% returns number of unique condition values of type conditionName
        function n = getnConditions(obj,conditionName)
            if ~ismember(conditionName,obj.conditions)
                disp([conditionName,' is not part of the conditions listed in the dataset']);
                n = 0;
                return
            end

            n = numel(unique(obj.fList.(conditionName)));
        end
        
        %% returns list of condition values of type conditionName
        function clist = getConditionsList(obj,conditionName)
            if ~ismember(conditionName,obj.conditions)
                disp([conditionName,' is not part of the conditions listed in the dataset']);
                clist = [];
                return
            end

            clist = unique(obj.fList.(conditionName));
        end
        
        %% converts select conditions to Numeric Format
        % conditions to be converted are listed in cell array conditionNames
        % if convertNonNumeralStrings = 0, it will not convert condition 
        % columns where some entries have non-numeric characters
        % if convertNonNumeralStrings = 1, it will remove all non-numeric
        % characters (everything except 0-9.) and convert the reminder to
        % numeric.
        function toNumeric(obj,conditionNames,convertNonNumeralStrings)
            c = intersect(conditionNames,obj.conditions);
            for i=1:numel(c)
                % check if the condition selected is convertible
                if ~isnumeric(obj.fList.(c{i}))
                    if ~iscell(obj.fList.(c{i}))
                        obj.fList.(c{i}) = cellstr(obj.fList.(c{i}));
                    end
                    
                    %check whether characters encode numbers only
                    clen = cellfun(@length,obj.fList.(c{i}),...
                        'UniformOutput',0);
                    cellfun(@num2str,clen,'UniformOutput',0);
                    regStr = cellfun(@strcat,...
                        repmat({'[0-9]{'},size(obj.fList,1),1),...
                        clen,...
                        repmat({'}'},size(obj.fList,1),1),...
                        'UniformOutput',false);
                    isNumeric = cellfun(@regexp,obj.fList.(c{i}),regStr,...
                        'UniformOutput',0);
                    isNumeric = sum(cellfun(@isempty,isNumeric)) == 0;

                    if isNumeric 
                        obj.fList.(c{i}) = cellfun(...
                            @str2double,obj.fList.(c{i}));
                    elseif convertNonNumeralStrings
                        % remove all non-numerals from expression and
                        % convert to string
                        regStr1 = repmat({'[^.0-9]'},size(obj.fList,1),1);
                        regStr2 = repmat({''},size(obj.fList,1),1);
                        obj.fList.(c{i}) = ...
                            cellfun(@regexprep,obj.fList.(c{i}),...
                            regStr1,regStr2,'UniformOutput',0);
                        obj.fList.(c{i}) = cellfun(...
                            @str2double,obj.fList.(c{i}));
                    end
                end
            end
        end
        
        %% converts select conditions from numeric to char Format
        % conditions to be converted are listed in cell array conditionNames
        % will not convert condition columns where some entries have
        % non-numeric characters. Can replace according to pattern: e.g. 
        % enter replacePattern = 'Round{}Hyb' to convert [1;2;3;] into
        % {'Round1Hyb';'Round2Hyb';'Round3Hyb';}
        % leave replacePattern empty otherwise.
        function toChar(obj,conditionNames,replacePattern)
            c = intersect(conditionNames,obj.conditions);
            for i=1:numel(c)
                if isnumeric(obj.fList.(c{i}))
                    obj.fList.(c{i}) = cellstr(num2str(obj.fList.(c{i})));
                    if ~isempty(replacePattern)
                        if contains(replacePattern,'{}')
                            k = strfind(replacePattern,'{}');
                            regStr1 = repmat({replacePattern(1:k-1)},...
                                size(obj.fList,1),1);
                            regStr2 = repmat({replacePattern(k+2:end)},...
                                size(obj.fList,1),1);
                            obj.fList.(c{i}) = cellfun(@strcat,...
                                regStr1,obj.fList.(c{i}),regStr2,...
                                'UniformOutput',false);
                        end
                    end
                end
            end
        end
        
        %% remove all file sets that are incomplete
        function removeIncompleteSets(obj)
            idxFileTypes = ismember(obj.fList.Properties.VariableNames,...
                obj.fileTypes);
            idxFileTypes = find(idxFileTypes);
            idxToKeep = [];
            for i=1:size(obj.fList,1)
                curFList = obj.getFileNameFromIdx(i,idxFileTypes);
                if ~ ismember(obj.voidEntryString, curFList)
                    idxToKeep = [idxToKeep;i];
                end
            end
            disp(['Removed ',num2str(size(obj.fList,1)-numel(idxToKeep)),'/',...
                num2str(size(obj.fList,1)),' rows with incomplete file sets.' ]);
            obj.fList = obj.fList(idxToKeep,:);
        end
                 
        %% add a filetype
        function addFileType(obj,fileTypeName,fileTypeEntries)
            if ismember(fileTypeName,obj.fileTypes)
                disp(['cannot add fileType ',fileTypeName,' to fList, ',...
                'it already exists']);
                return;
            end
            
            if isempty(fileTypeEntries)
                fileTypeEntries = repmat({obj.voidEntryString},...
                    size(obj.fList,1), 1);
            end
            
            if size(fileTypeEntries,1) ~= size(obj.fList,1) ...
                    || size(fileTypeEntries,2) ~= 1
                disp(['cannot add fileType ',fileTypeName,' to fList, ',...
                'entries have the wrong format']);
                return;
            end
            
            obj.fList.(fileTypeName) = fileTypeEntries;
            obj.fileTypes = [obj.fileTypes,fileTypeName];
        end
        
        %% delete a filetype 
        function deleteFileType(obj,fileTypeName)
            obj.fileTypes = setdiff(obj.fileTypes,fileTypeName);
            
            if ismember(fileTypeName,obj.fList.Properties.VariableNames)
                obj.fList = removevars(obj.fList,{fileTypeName});
            else
                disp(['Cannot remove file type ',fileTypeName,...
                    ' from fList, it is not present.']);
            end
            
        end
        
        %% check the value of the batchMode parameter in the config file
        function isbatch = checkBatchModeFromConfig(obj)
            
            % load ini config file
            inConf = obj.loadIni();
            
            % get batchMode status
            p = GetValues(inConf, ...
                'FileOptions', 'batchMode', 'fail');
            if strcmp(p,'fail') 
                isbatch = 0;
                return
            else
                isbatch = p;
                obj.options.batchMode = p;
            end
        end
        
        %% duplicate fileSet
        function fs = duplicate(obj)
            fs = fileSet;
            pNames = properties(fileSet);
            for i=1:numel(pNames)
                fs.(pNames{i}) = obj.(pNames{i});
            end
        end
        
        %% removes any rows in the file list that share the same set of conditions
        function cleanDuplicateConditions(obj)
            for i=size(obj.fList):-1:2
                curProperties = obj.fList(i,...
                    ismember(obj.fList.Properties.VariableNames,...
                    obj.conditions));
                prevProperties = obj.fList(1:i-1,...
                    ismember(obj.fList.Properties.VariableNames,...
                    obj.conditions));
                 if ismember(curProperties,prevProperties)
                     idx = setdiff(1:size(obj.fList),i);
                     obj.fList = obj.fList(idx,:);
                 end
            end
        end
        
        %% returns the file path of any file matching a desired type and 
        % set of conditions (output is a cell array of names)
        function [fname,idxConditionVals,idxFileType] = getFileName(obj,...
            conditionNames,conditionValues,desiredFileType,varargin)
        
            % obj: the fileSet object 
            % condition names: cell array listing the list of conditions names that are
            % constrained, e.g. {'FOV','Round'}
            % condition values: cell array lilsting the list of conditions values that
            % are constrained, they should match the condition names, eg. {1,2}
            % desiredFileType: the name of the variable encoding the location to the
            % desired file type, e.g. 'TransMat' for the column encoding transformation
            % Matrix

            % optional argument: 'first' or 'all' that decides whether to return the
            % first file filling the required criteria, or all of them (default) in
            % case there are more than one file that satisfy the conditions

            % example use: 
            % fname = getFileName(obj.fList,{'FOV','Round'},{2,3},'TransMat','all')
            % will return names of all Transformation Matrices that correspond to FOV 2 and round 3
            % as a cell array.

            if numel(conditionNames) ~= numel(conditionValues)
                fname = [];
                disp(['could not load file name from table - condition Names and ',...
                    'condition Values have different number of entries.']);
                return
            end

            % figure out whether or not to return all files filling the
            % condition
            expectedProcessMultiple = {'all','first'};
            defaultProcessMultiple = 'all';
            if numel(varargin)>0
                processMultiple = varargin{1};
                if ~ismember(processMultiple,expectedProcessMultiple)
                   processMultiple = defaultProcessMultiple; 
                end
            else
                processMultiple = defaultProcessMultiple; 
            end

            % locate position of condition variables among columns
            idxCon = zeros(numel(conditionNames));
            for i=1:numel(conditionNames)
                x = ismember(obj.conditions,conditionNames{i});
                if sum(x) == 0
                    disp(['could not load file name from table - required condition ',...
                        conditionNames{i},' is missing from condition list.']);
                    fname = [];
                    return;
                else
                    x = ismember(obj.fList.Properties.VariableNames,...
                        conditionNames{i});
                    idxCon(1,i) = find(x,1);
                end
            end

            % locate position of desired fileType
            x = ismember(obj.fileTypes,desiredFileType);
            if x == 0
                disp(['could not load file name from table - required file type ',...
                    desiredFileType,' is missing from fileTypes list.']);
                fname = [];
                return;
            else
                x = ismember(obj.fList.Properties.VariableNames,desiredFileType);
                idxFileType = find(x,1);
            end

            % find conditions that match the requested values
            idxConditionVals = true(size(obj.fList,1),1);
            for i=1:numel(conditionNames)
                if ischar(conditionValues{i})
                    conditionValues{i} = convertCharsToStrings(...
                        conditionValues{i});
                end
                t = table(conditionValues{i},'VariableNames',conditionNames(i));
                idxConditionVals = idxConditionVals & ismember(obj.fList(:,idxCon(1,i)), t );
            end

            if sum(idxConditionVals) == 0
                disp(['Could not find a ',desiredFileType,...
                    ' for requested round and FOV combination.']);
                fname = [];
                return;
            end

            % collect only the first file if option is selected, otherwise all files
            % will be returned
            if sum(idxConditionVals) > 1 && strcmp( processMultiple, 'first')
                idxConditionVals = find(idxConditionVals,1);   
            else
                idxConditionVals = find(idxConditionVals);
            end

            fname = obj.getFileNameFromIdx(idxConditionVals,idxFileType);
        end
   
        %% set file names matching a desired type and 
        % set of conditions 
        function setFileName(obj,...
            conditionNames,conditionValues,desiredFileType,fileNameEntry)
            % obj: the fileSet object 
            % condition names: cell array listing the conditions names that are
            % constrained, e.g. {'FOV','Round'}
            % condition values: cell array listing the conditions values that
            % are constrained, the size and order should match the condition names, eg. {1,2}
            % desiredFileType: the name of the desired file type, e.g. 'Loc' 
            % fielNameEntry: the filename to be inserted in the cells
            % matching the condition/fileType values.

            if numel(conditionNames) ~= numel(conditionValues)
                disp(['could not fill file name - condition Names and ',...
                    'condition Values have different number of entries.']);
                return
            end
            
            if ~ischar(fileNameEntry) && ~isstring(fileNameEntry)
                disp(['could not fill file name - file name should ',...
                    'be a string or char.']);
                return
            end
            
            % locate position of condition variables among columns
            idxCon = zeros(numel(conditionNames));
            for i=1:numel(conditionNames)
                x = ismember(obj.conditions,conditionNames{i});
                if sum(x) == 0
                    disp(['could not fill file name - required condition ',...
                        conditionNames{i},' is missing from condition list.']);
                    fname = [];
                    return;
                else
                    x = ismember(obj.fList.Properties.VariableNames,...
                        conditionNames{i});
                    idxCon(1,i) = find(x,1);
                end
            end

            % locate position of desired fileType
            x = ismember(obj.fileTypes,desiredFileType);
            if x == 0
                disp(['could not fill file name - required file type ',...
                    desiredFileType,' is missing from fileTypes list.']);
                return;
            else
                x = ismember(obj.fList.Properties.VariableNames,desiredFileType);
                idxFileType = find(x,1);
            end

            % find conditions that match the requested values
            idxConditionVals = true(size(obj.fList,1),1);
            for i=1:numel(conditionNames)
                if ischar(conditionValues{i})
                    conditionValues{i} = convertCharsToStrings(...
                        conditionValues{i});
                end
                t = table(conditionValues{i},'VariableNames',conditionNames(i));
                idxConditionVals = idxConditionVals & ismember(obj.fList(:,idxCon(1,i)), t );
            end

            if sum(idxConditionVals) == 0
                disp(['could not fill file name - Could not find a ',...
                    desiredVarName,...
                    ' for requested round and FOV combination.']);
                return;
            end
            idxConditionVals = find(idxConditionVals);
            obj.setFileNameFromIdx(idxConditionVals,idxFileType,fileNameEntry);
        end
              
        % gets file name(s) that match indices (output is a cell array)
        % (formalized as a method so it can be overriden in spotData where
        % the entries of the fList table are more complex).
        % returns a cell array
        function fname = getFileNameFromIdx(obj,idxConditionVals,idxFileType)
            if obj.areFListIndicesValid(idxConditionVals,idxFileType)
                fname = obj.fList{idxConditionVals,idxFileType};
            else
                disp('Cannot get File Name from fList, some indices out of bound');
                fname = {};
            end
        end
        
        % sets file name(s) that match indices 
        % (formalized as a method so it can be overriden in spotData where
        % the entries of the fList table are more complex).
        function setFileNameFromIdx(obj,idxConditionVals,idxFileType,fName)
            if obj.areFListIndicesValid(idxConditionVals,idxFileType)
                obj.fList{idxConditionVals,idxFileType} = {fName};
            else
                disp('Cannot set File Name in fList, some indices out of bound');
            end
        end
        
                %% initializes and populates obj options from config file
        function status = loadOptionsFromConfig(obj)
            
            % check whether ini file exists and load it
            inConf = obj.loadIni();
            if isempty(inConf)
                disp(['Could not load fileSet from config file: ',...
                    obj.iniFileName]);
                status = 0;
                return
            else 
                status = 1;
            end
            
            if isempty(obj.options)
                obj.initOptions;
            end
            
            fn = fieldnames(obj.options);
            for k=1:numel(fn)
                
                p = GetValues(inConf, ...
                    'FileOptions', fn{k}, 'fail');
                if strcmp(p,'fail') 
                    if isfield(obj.defaultOptions,fn{k})
                        obj.options.(fn{k}) = obj.defaultOptions.(fn{k});
                    else
                        status = 0;
                        obj.options.(fn{k}) = [];
                    end
                else
                     obj.options.(fn{k}) = p;
                end
            end
        end
        
        %% loads conditions from config
        function status = loadConditionsFromConfig(obj)
            
            % check whether ini file exists and load it
            inConf = obj.loadIni();
            if isempty(inConf)
                disp(['Could not load fileSet from config file: ',...
                    obj.iniFileName]);
                status = 0;
                return
            else
                status = 1;
            end
            
            % capture condition names from file name patterns
            [keys, count_keys] = inConf.GetKeys('Files');
            loadedConditions = {};
            for i=1:count_keys
                p = GetValues(inConf,'Files', keys{i}, 'fail');
                if ~strcmp(p,'fail') 
                    cd = obj.collectConditionsFromFilePattern(p);
                    loadedConditions = union(loadedConditions, cd);
                    loadedConditions = reshape(loadedConditions,...
                        1,numel(loadedConditions));
                end
            end
            
            if iscell(obj.conditions)
                obj.conditions = union(obj.conditions,loadedConditions);
                obj.conditions = reshape(obj.conditions,...
                    1,numel(obj.conditions));
            else
                obj.conditions = loadedConditions;
            end
            
            if ~obj.options.batchMode && isempty(obj.conditions)
                % set default condition name for single image
                obj.conditions = {obj.defaultConditionNameSingleImg}; 
            end
            
            % display list of conditions
            if ~isempty(obj.conditions)
                str = 'Loaded the following condition(s) from config file: ';
                for i=1:numel(obj.conditions)
                    str = [str,' ',obj.conditions{i},';'];
                end
                str(end) = '.';
                disp(str);
            end
        end
        
        %% loads file types from config
        function status = loadFileTypesFromConfig(obj)
            % check whether ini file exists and load it
            inConf = obj.loadIni();
            if isempty(inConf)
                disp(['Could not load fileSet from config file: ',...
                    obj.iniFileName]);
                status = 0;
                return
            else
                status = 1;
            end
            
            % parse list of file Types from ini config file
            obj.fileTypes = inConf.GetKeys('Files');
            obj.fileTypes = reshape(obj.fileTypes,1,numel(obj.fileTypes));
            
            if isempty(obj.fileTypes)
                disp(['Error: Could not parse the list of file Types in config file: ',...
                    obj.iniFileName,';']);
                disp(['    -> make sure the section [Files] is populated ',...
                    'in config file and contains file pattern entries',...
                    ' following format in example below:']);
                disp('    [Files]');
                disp('    loc = Path/to/Data/{Channel}/*_{Treatment}_{FOV}.loc3');
                disp('    mask = Path/to/Data/{Channel}/*_{Treatment}_{FOV}.tif');
                disp('    ...');
                return
            end
            % display list of fileTypes
            str = 'Loaded the following fileType(s) from config file: ';
            for i=1:numel(obj.fileTypes)
                str = [str,' ',obj.fileTypes{i},';'];
            end
            str(end) = '.';
            disp(str);
        end
        
        %% loads folder info from config
        function status = loadOutFolderInfoFromConfig(obj)
            
            % check whether ini file exists and load it
            inConf = obj.loadIni();
            if isempty(inConf)
                disp(['Could not load fileSet from config file: ',...
                    obj.iniFileName]);
                status = 0;
                return
            else
                status = 1;
            end
            
            % get output folder status
            p = GetValues(inConf, ...
                'Output', 'outFolder', 'fail');
            if ~strcmp(p,'fail') 
                obj.outFolder = p;
            else
                disp('Could not parse output folder from config file!');
            end
        end
        
        %% load patterns of file names from ini config. 
        % config filename is inputFileName (leave empty to use internally
        % stored obj.iniFileName)
        function status = loadPatternSetFromConfig(obj)
            
            % initialize patternSet
            obj.initPatternSet;
            
            % load file path patterns for each file type (localization, ...)
            for i=1:numel(obj.fileTypes)
                obj.loadSinglePatternFromConfig(obj.fileTypes{i});
            end
            
            % initialize condition list
            obj.initfList;
            
            status = 1;
        end
        
        %% loads the file path pattern for a specific file type into
        %fset.patternSet.([fileType,'_Filename'])
        function status = loadSinglePatternFromConfig(obj,fileType)
            
            % fileType can only be one of the allowed types
            if ~ismember(fileType,obj.fileTypes)
                disp(['Error: Cannot load fileType pattern; ',fileType,...
                    ' not found in fileTypes list.']);
                status = 0;
                return
            end 
            
            % check whether ini file exists and load it
            inConf = obj.loadIni();
            if isempty(inConf)
                disp(['Could not load fileSet from config file: ',...
                    obj.iniFileName]);
                status = 0;
                return
            else
                status = 1;
            end
            
            % making sure the fileType pattern is populated in the ini object
            p = GetValues(inConf,'Files', fileType, 'fail');
            if strcmp(p,'fail') 
                disp(['Error: Could not find ',fileType,...
                    ' naming pattern in config file: ',...
                    obj.iniFileName]);
                disp(['    -> make sure the section [Files] is populated ',...
                    'in config file and contains file pattern entry',...
                    ' following format in example below:']);
                disp('    [Files]');
                disp(['    ',fileType,' = Path/to/Data/{Channel}/*_{Treatment}_{FOV}.xyz']);
                disp('    ...');
                status = 0;
                return
            end
            
            
            % initialize patternSet if needed
            if isempty(obj.patternSet)
                obj.initPatternSet;
            end
            
            % add fileType to patternSet if needed
            if ~ismember(fileType,obj.patternSet.Properties.VariableNames)
                obj.patternSet = addvars(obj.patternSet,{0},...
                    'NewVariableNames',fileType);
            end
            
            % populate patternSet entry
            obj.patternSet.(fileType){1} = p;
        end
        
        %% assemble obj fList 
        function status = assembleFileListFromPattern(obj)
 
            % if not batch, populate the file list with the patternSet
            % entry (which should be a regular filename)
            if ~obj.checkBatchModeFromConfig
                obj.fillDefaultEntriesSingleFileMode;
                for i=1:numel(obj.fileTypes)
                    obj.fList.(obj.fileTypes{i}) = ...
                        obj.patternSet.(obj.fileTypes{i});
                end
                
            else
                % if batch mode, collect all files that match the pattern
                obj.fList = [];
                obj.conditions = [];
                for i=1:numel(obj.fileTypes)
                    [ftmp,ctmp] = ...
                        collectFilesFromSingleFileType(obj, obj.fileTypes{i});
                    [obj.fList, obj.conditions] = ...
                        obj.mergeFileLists(obj.fList,obj.conditions,ftmp,ctmp,1,1);
                end
            end
            status = 1;
        end
        
        %% merges two file lists
        function [fout,condNames] = mergeFileLists(obj,f1,c1,f2,c2,...
                includeIncompleteSets,mergeFileTypes)
        % merges the file lists contained in tables f1 and f2, based on their
        % respective sets of conditions c1 and c2. c1 and 2 should each be a subset
        % cell array of the variable names of resp. f1 and f2.

        % the script identifies the shared condition variables between f1 and f2,
        % then figures out the shared sets of conditions and matches entries in f1
        % and f2 by matching conditions.

        % includeIncompleteSets = 0/1: includes conditions with no match in the other file list.
        % empty entries are filled with standard missing tag: '<missing>' or
        % NaN for numeric types.
        
        % mergeFileTypes = 0/1: if the same fileType (= variable that is not a condition) 
        % is found in both f1 and f2, 
        % mergeFileTypes = 0 will create two separate file types 
            % in the output table with respective variable names 
            % <fileType>_1 and <fileType>_2;
        % mergeFileTypes = 1 will merge the two columns into one. In case
            % of shared conditions/filetypes combos, f1 entry is chosen.
        
        % if entries have no match or are missing from the input file tables, they
        % are populated with standard missing tags ('<missing>' or NaN depending on
        % whether the variable is numeric or not.

        % if uniqueMissingTag = 1, a unique numeric tag is added to each missing
        % entry. If = 0, all  missing entries share the same string.    
            
            % if one or more of the tables are empty
            if ~istable(f1)
                if obj.options.verbose
                    disp('file list f1 is not a table.');
                end
                if ~istable(f2)
                    if obj.options.verbose
                        disp('Error: file list f2 is not a table, cannot merge');
                    end
                    fout = [];
                    condNames = union(c1,c2,'stable');
                    return
                else
                    fout = f2;
                    c = f2.Properties.VariableNames;
                    condNames = union(c1,c2,'stable');
                    condNames = intersect(c,condNames,'stable');
                    condNames = reshape(condNames,1,numel(condNames));
                    return
                end
            end

            % if condition is empty, there should be only one file in the
            % list
            if isempty(c1)
                if size(f1,1) > 1
                    if obj.options.verbose
                        disp('condition list c1 is empty but table has mutliple rows, cannot merge.');
                    end
                    fout = [];
                    condNames = [];
                    return
                end
            end

            if isempty(c2)
                if size(f2,1) > 1
                    if obj.options.verbose
                        disp('condition list c2 is empty but table has mutliple rows, cannot merge.');
                    end
                    fout = [];
                    condNames = [];
                    return
                end
            end

            % making sure conditions listed are actually variableNames in their respective table
            c1 = intersect(f1.Properties.VariableNames,c1,'stable');
            c2 = intersect(f2.Properties.VariableNames,c2,'stable');

            % cleanup missing entries in condition columns and replace with
            % voidEntryString
            f1 = cleanup_table_columns(f1,c1,obj.voidEntryString,obj.options.uniqueMissingTag);
            f2 = cleanup_table_columns(f2,c2,obj.voidEntryString,obj.options.uniqueMissingTag);

            % merge variable list, initialize empty output table
            condNames = union(c1,c2,'stable');
            conds = cell2table(cell(0,numel(condNames)),'VariableNames',condNames);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            % populate all conditions
            % list of common variables
            cCommon = intersect(c1,c2);

            % extract part of table f1 that has variables common with f2
            [isCommon1,isUnique1,fCommon1,fUnique1] = ...
                obj.extractCommonConditions(f1,c1,cCommon);
            % same for f2
            [isCommon2,isUnique2,fCommon2,fUnique2] = ...
                obj.extractCommonConditions(f2,c2,cCommon);

            % build placeholders rows with variables compatible with f1 and f2 
            % to fill in for non-complete sets
            f1fill = obj.generatePlaceHolderTable(fUnique1);
            f2fill = obj.generatePlaceHolderTable(fUnique2);

            % find position index of each condition from list condNames in
            % the table comprising <common condition><unique conditions in f1><unique conditions in f2>
            order1 = zeros(1,numel(condNames));
            for j=1:numel(condNames)
                order1(j) = find(...
                    ismember( [ f1.Properties.VariableNames(isCommon1),...
                    f1.Properties.VariableNames(isUnique1),...
                    f2.Properties.VariableNames(isUnique2)],condNames{j} ));
            end

            % find position index of each condition from list condNames in
            % the table comprising <common condition><unique conditions in f2><unique conditions in f1>
            order2 = zeros(1,numel(condNames));
            for j=1:numel(condNames)
                order2(j) = find(...
                    ismember( [ f2.Properties.VariableNames(isCommon2),...
                    f2.Properties.VariableNames(isUnique2),...
                    f1.Properties.VariableNames(isUnique1)],condNames{j} ));
            end

            % match each entry of f1 with entries in f2 that share the same
            % set of common conditions
            n1 = size(f1,1);
            originalIdx = zeros(0,2); % will store indices of each file in the original table to track back file names later.
            for i=1:n1
                % find all members of f2 that share the same common
                % conditions with ith entry in f1
                idx2 = ismember(fCommon2,fCommon1(i,:));
                n2 = sum(idx2);

                % fill with place holder if needed
                if n2 == 0 && includeIncompleteSets
                    ftmp = [[fCommon1(i,:),fUnique1(i,:)],f2fill];
                    idx2 = NaN;
                else
                    ftmp = [repmat([fCommon1(i,:),fUnique1(i,:)],n2,1),...
                        fUnique2(idx2,:)];
                    idx2 = (find(idx2));
                end
                conds = [conds;ftmp(:,conds.Properties.VariableNames(order1))];

                originalIdx = [originalIdx; ...
                    [repmat(i,size(ftmp,1),1) , idx2  ] ];
            end

            % match each entry of f2 with entries in f1 that share the same
            % set of common conditions
            n2 = size(f2,1);
            for i=1:n2
                % find all members of f2 that share the same common
                % conditions with ith entry in f1
                idx1 = ismember(fCommon1,fCommon2(i,:));
                n1 = sum(idx1);

                if n1 == 0 && includeIncompleteSets
                    ftmp = [[fCommon2(i,:),fUnique2(i,:)],f1fill];
                    idx1 = NaN;
                else
                    ftmp = [repmat([fCommon2(i,:),fUnique2(i,:)],n1,1),fUnique1(idx1,:)];
                    idx1 = (find(idx1));
                end
                conds = [conds;ftmp(:,conds.Properties.VariableNames(order2))];
                originalIdx = [originalIdx; ...
                    [idx1 , repmat(i,size(ftmp,1),1) ] ];
            end

            % remove duplicate condition sets
            [conds, uniqueIdx] = unique(conds,'stable','rows');
            originalIdx = originalIdx(uniqueIdx,:);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % populate files corresponding to each set of conditions
            
            % extract the non-condition columns from each table
            f1out = f1(:,~ismember(f1.Properties.VariableNames, c1 ));
            f2out = f2(:,~ismember(f2.Properties.VariableNames, c2 ));
            
            % change names of variables to avoid duplicate names
            if ~mergeFileTypes
                repeatVars1 = intersect(f1out.Properties.VariableNames,...
                    [conds.Properties.VariableNames,f2out.Properties.VariableNames]);
                for i=1:numel(repeatVars1)
                    i1 = ismember( f1out.Properties.VariableNames,repeatVars1{i} );
                    f1out.Properties.VariableNames{ i1 } = ...
                        [ f1out.Properties.VariableNames{ i1 }, '_1'];
                end

                repeatVars2 = intersect(f2out.Properties.VariableNames,...
                    [conds.Properties.VariableNames,f1out.Properties.VariableNames]);
                for i=1:numel(repeatVars2)
                    i2 = ismember( f2out.Properties.VariableNames,repeatVars2{i} );
                    f2out.Properties.VariableNames{ i2 } = ...
                        [ f2out.Properties.VariableNames{ i2 }, '_2'];
                end
            end
            % split f1 and f2 into shared fileTypes and fileTypes unique to
            % each table
            sharedConditions = intersect(...
                f1out.Properties.VariableNames,...
                f2out.Properties.VariableNames);
            
            f1shared = f1out(: , ...
                ismember(f1out.Properties.VariableNames,sharedConditions));
            f2shared = f2out(: , ...
                ismember(f2out.Properties.VariableNames,sharedConditions));
            f1unique = f1out(: , ...
                ~ismember(f1out.Properties.VariableNames,sharedConditions));
            f2unique = f2out(: , ...
                ~ismember(f2out.Properties.VariableNames,sharedConditions));
            
            % generate place holder table rows for condition sets with missing entries
            f1unique_ph = obj.generatePlaceHolderTable(f1unique);
            f2unique_ph = obj.generatePlaceHolderTable(f2unique);

            fout = [];
            for i=1:size(originalIdx,1)
                if isnan(originalIdx(i,1))
                    f1unique_tmp = f1unique_ph;
                    f1shared_tmp = [];
                else
                    f1unique_tmp = f1unique(originalIdx(i,1),:);
                    f1shared_tmp = f1shared(originalIdx(i,1),:);
                end

                if isnan(originalIdx(i,2))
                    f2unique_tmp = f2unique_ph;
                    f2shared_tmp = [];
                else
                    f2unique_tmp = f2unique(originalIdx(i,2),:);
                    if isnan(originalIdx(i,1))
                        f2shared_tmp = f2shared(originalIdx(i,2),:);
                    else
                        f2shared_tmp = [];
                    end
                end

                fout = [fout;conds(i,:),f1shared_tmp,f2shared_tmp,...
                    f1unique_tmp,f2unique_tmp];
            end
            condNames = conds.Properties.VariableNames;
        end
        
    end
    
    %%
    methods (Access = 'protected')
        
        %% collect files in input folder that fit the pattern associated to 
        % one  fileType (entry for fileType must be populated in obj.patternSet)
        function [flist,conditions] = collectFilesFromSingleFileType(obj, fileType)
 
            % fileType can only be one of the allowed types
            if ~ismember(fileType,obj.fileTypes)
                disp(['Error: File type ',fileType,' is not allowed']);
                flist = [];
                conditions = [];
                return
            end 
            
            if ~ismember(fileType,...
                    obj.patternSet.Properties.VariableNames)
                disp(['Error: Pattern for ',fileType,' not found in table']);
                flist = [];
                conditions = [];
                return
            end
            
            % collect pattern for fileType from patternSet table 
            p = obj.patternSet.(fileType){1};
            
            [flist,conditions] = getFilesThatFollowPattern(p,fileType,...
                    obj.options.recursive,obj.options.useTerminal,...
                    obj.options.verbose);
                
        end
        
        %% load a list from the entry {iniSectionName,iniKeyName} in the ini config file, parses it as
        % a cell array and saves it into the object property
        % 'destinationPropertyName'     
        function status = readListEntryFromConfig(obj,destinationPropertyName,...
            iniSectionName,iniKeyName)
            removePaddingWhiteSpaces = 1;
            
            % load ini config file
            inConf = obj.loadIni();
            
            % get batchMode status
            p = GetValues(inConf, ...
                iniSectionName, iniKeyName, 'fail');
            if strcmp(p,'fail') 
                status = 0;
                return
            else
                 p = obj.parseCommaSeparatedList(p,...
                     removePaddingWhiteSpaces);
            end
            
            % if no entry is read, set to empty
            if numel(p) == 1 && isempty(p{1})
                p = [];
            end
            
            obj.(destinationPropertyName) = p;
            status = 1;
        end
        
        %% fill the conditions to default values when batch mode isnt
        % selected
        function fillDefaultEntriesSingleFileMode(obj)
            for i=1:numel(obj.conditions)
                if size(obj.fList,1) == 0
                    obj.fList = [obj.fList;repmat({''},1,size(obj.fList,2))];
                end
                obj.fList.(obj.conditions{i}){1} = ...
                        obj.defaultConditionValueSingleImg;
            end
        end
        
        %% initialize options and set to defaults
        function initOptions(obj)
            obj.options.recursive = obj.defaultOptions.recursive; 
            obj.options.verbose = obj.defaultOptions.verbose;
            obj.options.useTerminal = obj.defaultOptions.useTerminal;
            obj.options.batchMode = obj.defaultOptions.batchMode;
            obj.options.uniqueMissingTag = obj.defaultOptions.uniqueMissingTag;
        end
        
        %% initialize options to defaults
        function initDefaultOptions(obj)
            obj.defaultOptions.recursive = 0;
            obj.defaultOptions.verbose = 0;
            obj.defaultOptions.useTerminal = 1;
            obj.defaultOptions.batchMode = 1;
            obj.defaultOptions.uniqueMissingTag = 1;
        end
        
        %% initialize the patternSet table with the right variable number and names
        % based on fileTypes entry
        function initPatternSet(obj)
            if isempty(obj.patternSet)
                obj.patternSet = cell2table(cell(1,numel(obj.fileTypes)),...
                    'VariableNames',obj.fileTypes);
            end
        end
        
        %% initialize the patternSet table with the right variable number and names
        % based on fileTypes entry
        function initfList(obj)
            varNames = [obj.conditions,obj.fileTypes];
            obj.fList = cell2table(cell(0,numel(varNames)),...
                    'VariableNames',varNames);
        end
        
        function status = areFListIndicesValid(obj,...
                idxConditionVals,idxFileType)
            
            if isempty(obj.fList)
                status = 0;
                return
            end
            s = size(obj.fList);
            if min(idxConditionVals(:))<1 || min(idxFileType(:))<1 ...
                    || max(idxConditionVals(:))>s(1) ...
                    || max(idxFileType(:))> s(2)
                status = 0;
            else
                status = 1;
            end
            
        end
        
        %% parse a string that is comma separated, removing padding white
        % spaces (optional)
        function plist = parseCommaSeparatedList(~,strlist,removePaddingSpaces)
            plist = strsplit(strlist,',');
            if removePaddingSpaces
                for i=1:numel(plist)
                    remainsSpace = 1;
                    while remainsSpace
                        if startsWith(plist{i},' ')
                            plist{i} = plist{i}(2:end);
                        else
                            remainsSpace = 0;
                        end
                    end
                    
                    remainsSpace = 1;
                    while remainsSpace
                        if endsWith(plist{i},' ')
                            plist{i} = plist{i}(1:end-1);
                        else
                            remainsSpace = 0;
                        end
                    end
                end
            end
        end
        
        %% checking whether ini config file exists
        function status = checkIniFile(obj)
            % making sure the ini file exists
            if isempty(obj.iniFileName)
                disp(['Could not find ini config file : ',...
                        obj.iniFileName]);
                status = 0;
                return
            else
                if ~exist(obj.iniFileName,'File')
                    disp(['Could not find ini config file : ',...
                        obj.iniFileName]);
                    status = 0;
                    return
                end
            end
            status = 1;
        end
        
        %% load ini config file, returns empty if fails
        function inConf = loadIni(obj)
            
            % making sure the ini file exists
            if ~obj.checkIniFile
                inConf = [];
                return
            end
            
            % load ini config file
            inConf = IniConfig();
            inConf.ReadFile(obj.iniFileName);
        end
        
        % find condition names from string (should be enclosed in {})
        % formats them by removing spaces, and field size (e.g. '{ Channel
        % # 3}' becomes 'Channel'
        function cd = collectConditionsFromFilePattern(~,p)
            [~, cd, ~] = formatPatternStringForRegExpSearch(p);
        end
        
        %% extract part of table f (subfunction of mergeFileLists)
        % (with list of conditions conditionList) that shares common 
        % conditions with list commonConds
        function [isCommon,isUnique,fCommon,fUnique] = ...
                extractCommonConditions(~,f,conditionList,commonConds)
            isCommon = ismember(f.Properties.VariableNames,commonConds);
            fCommon = f(:, f.Properties.VariableNames(isCommon));
            isUnique =  ismember(f.Properties.VariableNames,...
                setdiff(conditionList,commonConds));
            fUnique = f(:, f.Properties.VariableNames(isUnique));
        end
        
        %% generate place holder table with same variables as input table f.
        % (subfunction of mergeFileLists) 
        % one row is filled with either standardized empty String or value
        % depending on the variable type.
        function ph = generatePlaceHolderTable(obj,f)
            if isempty(f)
                ph = [];
                return
            end
            for i=1:size(f,2)
                if isnumeric(f.(f.Properties.VariableNames{i}))
                    if i==1
                        ph = table(obj.voidEntryValue,...
                            'VariableNames',f.Properties.VariableNames(i));
                    else
                        ph = addvars(ph,obj.voidEntryValue,...
                            'NewVariableNames',f.Properties.VariableNames(i));
                    end
                else
                    if i==1
                        ph = table({obj.voidEntryString},...
                            'VariableNames',f.Properties.VariableNames(i));
                    else
                        ph = addvars(ph,{obj.voidEntryString},...
                            'NewVariableNames',f.Properties.VariableNames(i));
                    end
                end
            end
        end
        
        
    end
end