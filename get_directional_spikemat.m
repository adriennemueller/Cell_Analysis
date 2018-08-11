% Helper function that returns a structured set of spikemats for different
% directions
function spikemat_struct = get_directional_spikemat( spikemat, current, window, inout, contrast_flag  )
    spikemat_struct = struct;
    directions = unique([spikemat.theta]);
    % Loop through directions and place into struct
    for i = 1:length(directions)
        spikemat_struct(i).direction = directions(i);
        spikemat_struct(i).spikes = filtered_windowed_spikemat( spikemat, current, window, directions(i), inout );
        
        % Add contrast data, if attend_Contrast paradigm
        if contrast_flag
            spikemat_struct(i).contrast = get_contrasts( spikemat, current, directions(i), inout );
        end
        
    end
        
end

