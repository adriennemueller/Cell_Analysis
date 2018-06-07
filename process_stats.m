% Will run through all processed files and calculate d's, anovas, etc for
% different epochs of the trial and save them to master_file_struct
function mfs = process_stats( mfs )

    if nargin < 1, load master_file_struct; mfs = master_file_struct; end
    

    % Loop through all processed files
    
    % Identify the paradigms present in those files
    
    % Identify whether there are enough correct trials within a given
    % paradigm to process it for stats
    
    % Calculate stats for each paradigm (anova, d', etc);
    
    % Save substruct of stats for each paradigm into mfs
    
    
    % Save out master_file_struct

end
