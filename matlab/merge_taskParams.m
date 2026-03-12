function merged_struct = merge_taskParams(tbl)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
all_params = {tbl.params.block}.';
all_params = cellfun(@convertSparse, all_params, 'uni', 0);

structStrings = cellfun(@(x) jsonencode(x), all_params, 'uni', 0);
[~, uniqueIdx] = unique(structStrings, 'stable');
unique_structs = all_params(uniqueIdx);
merged_struct = struct();
fieldNames = fieldnames(unique_structs{1});
for ii = 1:numel(fieldNames)
    field = fieldNames{ii};
    merged_struct.(field) = cellfun(@(s) s.(field), unique_structs, 'uni', 0);
    if all(cellfun(@isnumeric, merged_struct.(field)))
        merged_struct.(field) = cell2mat(merged_struct.(field));
    end
end
end

function s = convertSparse(s)
fields = fieldnames(s);

for k = 1:numel(fields)
    v = s.(fields{k});
    
    if issparse(v)
        s.(fields{k}) = full(v);
    end
    
end
end