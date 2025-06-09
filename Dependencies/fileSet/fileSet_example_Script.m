% example uses

% dependencies:
% ini config
% Evgeny Pr (2020). INI Config 
% (https://www.mathworks.com/matlabcentral/fileexchange/24992-ini-config), 
% MATLAB Central File Exchange. Retrieved May 22, 2020.

%% generate the dummy data set 
% also creates the variable config_fname that stores the path to the dataset config file
generate_example_dataSet;

%% building a fileSet object
%intialize object
fs = fileSet;

% populate object properties from ini config file
fs.buildFileSetFromConfig(config_fname);

% fs.fList lists all files in a table where the columns list each
% condition (Condition | Well | Channel | Time ), and each file type ( Img | WellMask).
% Each row lists a set of conditions and the files of each type that correspond to it.

%% example: manipulating a file set table

% convert entries in conditions columns from characters to numbers.

% the convertNonNumeralStrings flag ensures that all columns including those 
% that mix letters and numbers are converted - e.g. the time points entries in the 
% table {'t1'; 't2'; 't3'; ...}  to [1;2;3;...]. 
disp('Converting Time and Channel Conditions to numeric format...');
convertNonNumeralStrings = 1; 
fs.toNumeric({'Time','Channel'},convertNonNumeralStrings);
disp('Done converting Time and Channel Conditions to numeric format...');
disp(' ');

% this reverses the operation on the time column, adding the leading 't'
% back.
disp('Converting Time back to character, with format t1, t2, ...');
fs.toChar({'Time'},'t{}');
disp(' ');

% this reverses the operation on the Channel column, without adding any
% letters.
disp('Converting Channel back to character...');
fs.toChar({'Channel'},[]);
disp(' ');

% this removes any rows that share the same set of condition entries
% (keeps the first entry)
disp('Clearing duplicate conditions...');
fs.cleanDuplicateConditions();
disp(' ');

% remove all rows of file list that have missing entiries
disp('Removing incomplete file sets...');
fs.removeIncompleteSets();
disp(' ');

% collect names of all Img files matching t = 3 and channel = 1 conditions
disp('Collecting all Img files from time point t3, Channel 1 in variable fnames1...');
fnames1 = fs.getFileName({'Time','Channel'},{'t3','1'},'Img','all');
disp(' ');

% same, but colecting only the first file only
disp('Collecting the first Img file from time point t3, Channel 1 in variable fnames2...');
fnames2 = fs.getFileName({'Time','Channel'},{'t3','1'},'Img','first');
disp(' ');

% add a new file type to the fileSet called someNewFileType 
% all rows are populated with the default empty entry ('<missing>'):
disp('Adding a new fileType called someNewFileType to fs...');
fs.addFileType('someNewFileType',[]);
disp(' ');

%%
% fill all entries in FileType Loc column that match the set of conditions
% here Time =3 and Channel = 1 with the filename newFileName
disp('fill all Img file entries for Time t3, Channel 1 with file name newFileName');
fs.setFileName({'Time','Channel'},{'t3','1'},'Img','newFileName');
disp(' ');

% self explanatory
disp('deleting fileType someNewFileType');
fs.deleteFileType('someNewFileType');

%% example: combining file sets that share same condition names and file type names
% duplicate fileSet; keep just one row and modify the condition
disp('Duplicating fs into fs2, adding a new condition to fs2...');
fs2 = fs.duplicate();
fs2.fList =fs2.fList(1,:);
fs2.fList.Condition{1} = 'some new condition';
disp(' ');

% combine fs2 and fs
disp('recombining fs and fs2...');
includeIncompleteSets = 1; % this only matters when combining file lists with different sets of conditions columns
mergeFileTypes = 1; % this ensures that the columns for each file Type are merged
fs.fList = fs.mergeFileLists(fs.fList,fs.conditions,fs2.fList,fs2.conditions,...
                includeIncompleteSets,mergeFileTypes);
disp(' ');

%% write file list to file 
%empty argument to use default fine name and location, otherwise enter
%desired path as argument: fs.saveFileList('path/to/output/file/table.txt');
disp('writing file list to file...');
fs.saveFileList([]);

%% generate fileSet object fs3 from saved table file
% needs to be a tab delimited file, with variable names as column headers

fs3 = fileSet();
conditions = ''; % you can specific with table columns are the conditions;
% if left empty, a command line prompt will request them.

outFolder = ''; % you can specify the output folder for this new object, or leave empty
pathToSavedTable = [fs.outFolder,filesep,'fileSetList.txt'];
fs3.buildFileSetFromTableInFile(pathToSavedTable,conditions,outFolder);