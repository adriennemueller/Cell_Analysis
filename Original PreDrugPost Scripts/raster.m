function rslt = raster( times, startTime, endTime )
    startTime = startTime * 1000;
    endTime = endTime * 1000;

    % Make empty vector for times
    dur = endTime - startTime;
    vec = zeros(1,dur);

    % Convert to ms
    times = floor(times * 1000);
    
    % Realign to '0';
    times = times - startTime;
    
    for i = 1:length(times) 
        vec(times(i)) = 1;
    end
    
    rslt = vec;
end