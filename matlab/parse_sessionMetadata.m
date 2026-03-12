function [metadata, eye_chan_labels, diode_chan_label, pupil_chan_label] = parse_sessionMetadata(session_path)

if isfile(fullfile(session_path,'metadata.json'))
    metadata = loadMetadataJSON(fullfile(session_path,'metadata.json'));

    eye_chan_labels = metadata.HEeye_VEeye_diode_pupil(1:2)';
    diode_chan_label = metadata.HEeye_VEeye_diode_pupil{3};
    pupil_chan_label = metadata.HEeye_VEeye_diode_pupil{4};

else
    [~,fn,~] = fileparts(session_path);
    metadata.sess_name = fn;
    metadata.probe_type = 'behavior';

    eye_chan_labels = {'10241','10242'};
    diode_chan_label = '10243';
    pupil_chan_label = '10244';
end

end
