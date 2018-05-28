function rslt = update_spikemat( spikemat, eventcodes, newSpikeVal )

    [first_idx, second_idx] = translate_eventcodes( eventcodes );

    %%% Make this update instead of overwrite
    spikemat(first_idx, second_idx) = newSpikeVal;
    % Consider giving a border of empty values so can see 
    surf(spikemat,'EdgeColor','None', 'facecolor', 'interp'); view(2);
    rslt = spikemat;

end


function [first_idx, second_idx] = translate_eventcodes( eventcodes )

    coords = [-7 -6 -5 -4 -3 -2 -1 1 2 3 4 5 6 7];
    
    x_val = map_ecode( eventcodes(1) );
    y_val = map_ecode( eventcodes(2) );

    % x = second index because it is the 'column' index
    % y = first index because it is the 'row' index
    second_idx = find( x_val == coords );
    first_idx  = find( y_val == coords );
    
end


% This function takes the weird monkeylogic values that ML currently sends,
% and that Plexon makes negative, and returns their actual value.
function val = map_ecode( ecode )

    val = ecode;

    %%% TMP %%%
%     if ecode == 255
%         val = 9;
%     end
end