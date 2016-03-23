% psp_kernel generates a kernel of length (t) (in ms) shaped like a postsynaptic 
% potential.
% tau_g is the time constant for the growth phase
% tau_d is the time constant for the decay phase
% The formula it uses is from Thompson et al (Schall senior) (1996):
% R(t) = (1 - exp(-t/tau_d))*exp(-t/tau_d). They chose tau_g and tau_d 
% based on physiological data of excitatory synapses from Sayer et al (1990).%
% These are therefore our chosen default values: tau_g = 1 ms and tau_d = 20 ms.

function Rt = psp_kernel( t, tau_g, tau_d )

    if nargin < 2
        tau_g = 1; %Default growth time constant: 1ms
    end
    if nargin < 3
        tau_d = 20; %Default decay time constant: 20ms
    end

    Rt = (1 - exp( -t ./ tau_g )) .* exp( -t ./ tau_d );
    
    % Scale so total of values of kernel == 1.
    Rt = Rt ./ sum(Rt);

end
