function [flist,conditions, len] = ...
    getFilesThatFollowPattern(pattern,fileType,recursive,useTerminal,verbose)

% dirName is the folder whether the function searches for files following the
    % pattern, e.g. '/path/to/my/directory/'
% pattern is a char that specifies the naming pattern, e.g. 
    % pattern =
    % '/path/to/my/directory/*/{Channel}someCommonStringInFileNames{Time #3}{FOV #4}.loc3'
    % use curly braces with a key inside (e.g. {Channel} or {Time}) to
    % indicate variable portions of the names.
    % use #n to indicate that a field is exactly n characters long.
    % * work as wildcards. 

% fileType (optional, leave empty if not needed): the name of the variable 
    % storing the file names, e.g. 'Loc'. if left empty, variable is name 'fileName'.  
% useTerminal: enter 1 to try to use a terminal awk script, 0 to use the
    % matlab-only commands
    % Terminal/awk is a faster way to search folders with large arborescence; if awk isn't
    % installed, the function will revert to the default (slower) matlab dir function.
% recursive: enter 1 to search all subfolders recursively, 0 for just the
    % root level.

%output flist is a table with variables as follows:
    % <key 1> ... <key n> <fileType>
    % each row stores the variable segment associated to each key as a
    % char and the full file name under fileType. e.g. in example above:
    %  Channel | Time | FOV | fileName
    %   'C1'     't0'   #1    '/path/to/my/directory/x/y/z/C1someCommonStringInFileNamest0#1.loc3'
    %   'C2'     't0'   #1    '/path/to/my/directory/x/y/z/C2someCommonStringInFileNamest0#1.loc3'
    %   ...

% conditions is a 1 x nkeys cell array that lists all the keys used

% len is a 1 x nkeys cell array that lists the sizes of the fields
% corresponding to each key, e.g. '*' for arbirtrary size, '{m}' for
% exactly m characters, '{mn}' for any size between m and n.

%% extract the lowest directory level that does not have a wild card
% find the position of the last occurence of a {condition} in the pattern
if ~contains(pattern,'{')
    k1 = length(pattern);
else
    k1 = strfind(pattern,'{');
    k1 = k1(1)-1;
end

% find the position of the last occurence of a wild card in the pattern
if ~contains(pattern,'*')
    k2 = length(pattern);
else
    k2 = strfind(pattern,'*');
    k2 = k2(1)-1;
end

% find the last directory separator (/ or \) before * or { and collect the 
% directory name as dirName where to look for files.
k = min([k1,k2]);
dirName = pattern(1:k);
if ~contains(pattern,filesep)
    disp(['Could not parse parent directory from file name pattern ',...
        pattern,...
        '; make sure the full path is included in the pattern.']);
else
    k = strfind(dirName,filesep);
    k = k(end);
    dirName = dirName(1:k);
end


%% get list of all files in input folder
flist = getFileListFromDir(dirName,recursive,useTerminal);
if verbose
    disp(['Exploring ',num2str(numel(flist)),...
        ' total files in input directory ',dirName,'...']);
end

%% translate the pattern into a regexp query expression
[sp, keys, len] = formatPatternStringForRegExpSearch(pattern);

if verbose
    dispStr = [];
    for i=1:numel(keys)
        dispStr = [dispStr,' ',keys{i},';'];
    end
    disp(['Parsed ',num2str(numel(keys)),' keys in pattern:',dispStr]);
    
end

%% search for files that match the pattern in the full list
tokens = cell(0,1+numel(keys));
nt = 0;
for i=1:numel(flist)
    t = regexp(flist{i},sp,'tokens');
    if ~isempty(t)
        nt = nt+1;
        tokens(nt,1:numel(keys)) = t{1};
        tokens{nt,end} = flist{i};
    end
end
if verbose
    disp(['Found ',num2str(nt),' files matching pattern.']);
end

%% convert to table format
if isempty(fileType)
    fileType = 'FileName';
end
flist = cell2table(tokens,'VariableNames',[keys,{fileType}]);
conditions = keys;

end