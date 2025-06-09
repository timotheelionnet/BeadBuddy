function [sp, keys, len] = formatPatternStringForRegExpSearch(pattern)

% takes as input a search pattern string, e.g. 
    % pattern =
    % '/path/to/my/directory/*/{Channel}someCommonStringInFileNames{Time #3}{FOV #2#4}.loc3'
    % use curly braces with a key inside (e.g. {Channel} or {Time}) to
    % indicate variable portions of the names.
    % use #n at the end of a field to indicate that it is exactly n characters long.
    % use #n#m at the end of a field to indicate that it is between n and m characters long.
    % * work as wildcards. 

    % returns the string to be used as a regexp query, e.g. from example above:
    % '/path/to/my/directory/.*/.*/(.*)someCommonStringInFileNames(.*)(.*)\.loc3'
    % keys are replaced by tokens to capture their individual expressions
    % via regexp.
    
    % keys is a cell array that stores the keys in the curly braces, e.g.
    % in the example above
    % keys = {'Channel','Time','FOV'} (leading or trailing spaces are removed)
    % len is a cell array that stores the lengths (encoded as a regexp) of
    % each field.

%% reformat pattern string, save keys

% find portions of the pattern in between pairs of curly braces
[s,e] = regexp(pattern,'{[^{}]*}');

% collect keys (text between curly braces)
keys = [];
rawkeys = [];
len = [];
for nk=1:numel(s)
     kstr = pattern(s(nk)+1:e(nk)-1);
     rawkeys{nk} = kstr;
     hashtagIdx = strfind(kstr,'#');
     if isempty(hashtagIdx)
         keys{nk} = kstr;
         len{nk} = '*';
     else
         if numel(hashtagIdx) == 1
            keys{nk} = kstr(1:hashtagIdx(1)-1);
            curlen = round(str2double(kstr(hashtagIdx(1)+1:end)));
            if curlen > 0
                len{nk} = ['{',num2str(curlen),'}'];
            else
                disp(['Could not parse length of field ',keys{nk} ]);
                len{nk} = '*';
            end
         elseif numel(hashtagIdx) >= 2
            keys{nk} = kstr(1:hashtagIdx(1)-1);
            curlen = round(str2double(kstr(hashtagIdx(1)+1:hashtagIdx(2)-1)));
            if curlen > 0
                len{nk} = ['{',num2str(curlen)];
                curlen = round(str2double(kstr(hashtagIdx(2)+1:end)));
                if curlen > 0
                    len{nk} = [len{nk},',',num2str(curlen),'}'];
                else
                    disp(['Could not parse length of field ',keys{nk} ]);
                    len{nk} = '*';
                end
            else
                disp(['Could not parse length of field ',keys{nk} ]);
                len{nk} = '*';
            end 
         end
     end
end

% add escape characters in front of any special characters
sp = regexptranslate('escape',pattern);

% add back any wild cards in the name
sp = strrep(sp,'\*','.*');

% replace {key} with place holder regular expression in the pattern
for nk=1:numel(s)
    srep = ['{',rawkeys{nk},'}'];
    srep = regexptranslate('escape',srep);
    sp = strrep(sp,srep,['(.',len{nk},')']);
end

% select out the part of the path that is one level up and above from the pattern
if ~startsWith(sp,regexptranslate('escape',filesep))
    sp = ['.*',filesep,sp];
end

% remove leading/trailing spaces in keys
for nk = 1:numel(keys)
    while startsWith(keys{nk},' ')
        keys{nk} = keys{nk}(2:end);
    end
    while endsWith(keys{nk},' ')
        keys{nk} = keys{nk}(1:end-1);
    end
end

% replace spaces by underscores
for nk = 1:numel(keys)
    keys{nk} = strrep(keys{nk},' ','_');
end

end