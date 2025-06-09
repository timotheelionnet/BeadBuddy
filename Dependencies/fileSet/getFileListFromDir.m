function flist = getFileListFromDir(dirName,recursive,useTerminal)

% dirName is the folder whether the function searches for files following the
    % pattern.
    
% this awk line is used if the terminal option is selected

%awkScript = '| awk ''/:$/&&f{s=$0;f=0}
%/:$/&&!f{sub(/:$/,"");s=$0;f=1;next} NF&&f{ print s"/"$0 }'''; %older
%syntax that does not get the root directory files correctly

if ~endsWith(dirName,filesep)
    dirName = [dirName,filesep];
end

awkScript = ['| awk ''!/:$/&&!f{s="',dirName,'";f=1;} /:$/&&f{s=$0;f=0} /:$/&&!f{sub(/:$/,"");s=$0;f=1;next} NF&&f{ print s"/"$0 }'''];

%% get list of files inside directory 
flist_found = 0;
if useTerminal
    % check whether awk is installed
    [ ~, x] = system('awk - XXX');
    if startsWith(x,'awk:')
        if recursive
            [~,flist] = system(['ls -p -R ''',dirName,''' ',awkScript]);
        else
            [~,flist] = system(['ls -p ''',dirName,''' ',awkScript]);
        end
        flist = strrep(flist,[filesep,filesep],filesep);
        flist = strsplit(flist,'\n');
        flist = flist';
        flist_found = 1;
    end
end

if ~flist_found
    if recursive
        flist = dir(fullfile(dirName,'**/*'));
    else
        flist = dir(dirName);
    end
    flist = flist(~[flist.isdir])';
    flist = cellfun(@fullfile, {flist.folder},{flist.name},'UniformOutput',0)';
end

