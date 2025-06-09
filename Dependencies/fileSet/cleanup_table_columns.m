function t = cleanup_table_columns(t,varNames,voidEntryString,uniqueMissingTag)

% cleans up specific columns of a table t.
% columns to be cleaned up are specified by their variable names varNames
% empty entries in cell array columns are replaced by a standard string
% voidEmptyString.
% if uniqueMissingTag = 1, a unique numeric tag is added to each missing
% entry. If = 0, all  missing entries share the same string.


for i=1:numel(varNames)
    if ~isnumeric(t.(varNames{i}))
        % replace empty arrays with standard string
        t.(varNames{i})(cellfun(@isempty,t.(varNames{i}))) = ...
                {voidEntryString};
            
        % replace stray numerical values with their equivalent
        % string
        t.(varNames{i})(cellfun(@isnumeric,t.(varNames{i}))) = ...
            cellfun(@num2str,...
            t.(varNames{i})(cellfun(@isnumeric,t.(varNames{i}))),...
            'UniformOutput',0);
        
        % add unique tag to missing entries if needed    
        if uniqueMissingTag    
            nEmpty = sum(ismember(t.(varNames{i}),{voidEntryString}  ));
            if nEmpty ~=0
                t.(varNames{i})(ismember(t.(varNames{i}),{voidEntryString}  )) = ...
                    cellfun(@strcat,repmat({voidEntryString},nEmpty,1),...
                    cellstr(cellfun(@num2str, num2cell((1:nEmpty)'))),'UniformOutput',0);
            end
                
        end
    end
end
