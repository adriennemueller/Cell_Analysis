% Convolves a spike train with a PSP kernel
function rslt = spike_desity_psp( spikes )

    % Build a 200ms kernel with the default tau_g and tau_d
    kernel = psp_kernel( 1:200 );

    rslt = conv( spikes, kernel, 'same' );

end
