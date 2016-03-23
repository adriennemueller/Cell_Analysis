% Convolves a spike train with a Gaussian of the given standard deviation.
function rslt = spike_density( spikes, StDevInMS )

    if nargin < 2
        StDevInMS = 100; %Default standard deviation = 100ms;
    end
    
    SamplesPerMS = 1;
    NumberOfStDevs = 3;
    WindowSize = StDevInMS * SamplesPerMS * NumberOfStDevs;
    
    G = gausswin( WindowSize, NumberOfStDevs );
    G = G / sum(G); % Rescale gaussian to integrate to 1
    G = G * SamplesPerMS * 1000; % Rescale to spikes/sec
    
    rslt = conv( spikes, G, 'same' );
end