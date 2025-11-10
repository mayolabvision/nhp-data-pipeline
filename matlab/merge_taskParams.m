function merged_struct = merge_taskParams(tbl)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
all_params = {tbl.params.block}.';
structStrings = cellfun(@(x) jsonencode(x), all_params, 'UniformOutput', false);
[~, uniqueIdx] = unique(structStrings, 'stable');
unique_structs = all_params(uniqueIdx);
merged_struct = struct();
fieldNames = fieldnames(unique_structs{1});
for ii = 1:numel(fieldNames)
    field = fieldNames{ii};
    merged_struct.(field) = cellfun(@(s) s.(field), unique_structs, 'UniformOutput', false);
    if all(cellfun(@isnumeric, merged_struct.(field)))
        merged_struct.(field) = cell2mat(merged_struct.(field));
    end
end
end