% spike_overlay overlays spike traces from the the matrix of spike times,
% and voltage values output by plexon. After loading the plexon .mat file,
% pass the matrix as input into spike_overlay.
% It returns both a figure handle: fprop, and an axes handle: aprop, so the
% plot can be integrated into other plot-views.

function [fprop, aprop] = spike_overlay( spike_mat )
    [m, n] = size(spike_mat);
    
    fprop = figure();
    hold on;
   
    for i = 1:m
    plot(spike_mat(i,2:n), 'k'); % The first column is a list of times of spikes; 2-n are recorded voltages for the spikes.      
    end

    hold off;
    
    xlabel( 'Time (AU)' );
    ylabel( 'Voltage (uV)' );
    
    aprop = gca;
    
end