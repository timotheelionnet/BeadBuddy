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
convertNonNumeralStrings = 1; 
fs.toNumeric({'Time','Channel'},convertNonNumeralStrings);

% this reverses the operation on the time column, adding the leading 't'
% back.
fs.toChar({'Time'},'t{}');

% this reverses the operation on the Channel column, without adding any
% letters.
fs.toChar({'Channel'},[]);

% this removes any rows that share the same set of condition entries
% (keeps the first entry)
fs.cleanDuplicateConditions();

% remove all rows of file list that have missing entiries
fs.removeIncompleteSets();

% collect names of all Img files matching t = 3 and channel = 1 conditions
fnamesToProcess1 = fs.collectFileName({'Time','Channel'},{3,1},'Loc','all');

% same, but colecting only the first file only
fnamesToProcess2 = fs.collectFileName({'Time','Channel'},{3,1},'Loc','first');

% add a new file type to the fileSet called someNewFileType 
% all rows are populated with the default empty entry ('<missing>'):
fs.addFileType('someNewFileType',[]);

% fill all entries in FileType Loc column that match the set of conditions
% here Time =3 and Channel = 1 with the filename newFileName
fs.fillFileName({'Time','Channel'},{3,1},'Loc','newFileName');

% self explanatory
fs.deleteFileType('someNewFileType');

%% example: combining file sets that share same condition names and file type names
% duplicate fileSet; keep just one row and modify the condition
fs2 = fs.duplicate();
fs2.fList =fs2.fList(1,:);
fs2.fList.Condition{1} = 'some new condition';

% combine fs2 and fs
includeIncompleteSets = 1; % this only matters when combining file lists with different sets of conditions columns
mergeFileTypes = 1; % this ensures that the columns for each file Type are merged
fs.fList = fs.mergeFileLists(fs.fList,fs.conditions,fs2.fList,fs2.conditions,...
                includeIncompleteSets,mergeFileTypes);
                        
%% write file list to file 
%empty argument to use default fine name and location, otherwise enter
%desired path as argument: fs.saveFileList('path/to/output/file/table.txt');
fs.saveFileList([]);