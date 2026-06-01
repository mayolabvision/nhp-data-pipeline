function short_id = shorten_sorter_path(sorter_path, key_dir)
% SHORTEN_SORTER_PATH  Hash a long SORTER_PATH string to a short, repeatable ID.
%
%   short_id = shorten_sorter_path(sorter_path)
%   short_id = shorten_sorter_path(sorter_path, key_dir)
%
%   sorter_path : the full SORTER_PATH string (e.g. the chain of hashes
%                 separated by dashes that SpikeInterface generates)
%   key_dir     : directory in which to write/update 'sorter_path_key.json'
%                 Defaults to the current working directory.
%
%   short_id    : first 16 hex characters of the SHA-256 of sorter_path.
%                 Guaranteed to be the same every time for the same input.
%
%   The mapping  short_id -> sorter_path  is appended to
%   <key_dir>/sorter_path_key.json so you can always recover the original.

    if nargin < 2 || isempty(key_dir)
        key_dir = pwd;
    end

    % --- compute SHA-256 and take first 16 hex chars ---
    import java.security.MessageDigest
    import java.math.BigInteger

    md  = MessageDigest.getInstance('SHA-256');
    raw = md.digest(uint8(sorter_path));
    % convert signed bytes to unsigned hex
    hex_str = sprintf('%02x', mod(double(raw), 256));
    short_id = hex_str(1:16);

    % --- update the key file ---
    key_file = fullfile(key_dir, 'sorter_path_key.json');

    if isfile(key_file)
        key_map = jsondecode(fileread(key_file));
    else
        key_map = struct();
    end

    % struct field names cannot start with a digit, so prefix with 'id_'
    field = ['id_', short_id];
    if ~isfield(key_map, field)
        key_map.(field) = sorter_path;
        fid = fopen(key_file, 'w');
        fprintf(fid, '%s', jsonencode(key_map, 'PrettyPrint', true));
        fclose(fid);
    end
end
