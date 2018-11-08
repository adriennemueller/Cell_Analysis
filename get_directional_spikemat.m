% Helper function that returns a structured set of spikemats for different
% directions
function spikemat_struct = get_directional_spikemat( spikemat, current, window_str, inout, contrast_flag )
    spikemat_struct = struct;
    directions = unique([spikemat.theta]);
    % Loop through directions and place into struct
    for i = 1:length(directions)
        spikemat_struct(i).direction = directions(i);
        [spikemat_struct(i).spikes spikemat_struct(i).contrasts] = filtered_windowed_spikemat( spikemat, current, window_str, directions(i), inout, contrast_flag );
    end
        
end

