% add_session will add the different sources of data from a recording
% session into the master-struct which contains the info for all recording
% sessions.

% sub_direc   = subdirectory containing the files for this session
% event_file  = STRING of a plexon eventfile (.mat)
% bhv_files = cell array of STRINGS of monkeylogic  behav_files (ending in .bhv)
% unit_files  = cell array of STRINGS plexon unit files (.mat). 

% tmp_struct = an optional struct to use if you don't want to load and
% append to the master_struct
function master_file_struct = add_session( sub_direc, event_file, bhv_files, unit_files, drug, tmp_struct )

    
    % Load Master Struct if a temporary struct isn't defined.
    if nargin < 6
       load( 'master_file_struct.mat', 'master_file_struct' );
    else
       master_file_struct = tmp_struct;
    end
    
    % Append Information
    if isfield(master_file_struct, 'session')
        session_count = length( master_file_struct.session );
    else
        master_file_struct.session.sub_direc = struct;
        session_count = 0;
    end
    
    % Test whether this directory has already been added to the
    % master_struct; if so: skip it.
    if max(strcmp( {master_file_struct.session.sub_direc}, sub_direc ))
        return
    end
    
    master_file_struct.session( session_count + 1 ).sub_direc    = sub_direc;
    master_file_struct.session( session_count + 1 ).event_file   = event_file;
    master_file_struct.session( session_count + 1 ).bhv_files    = bhv_files;
    master_file_struct.session( session_count + 1 ).unit_files   = unit_files;
    master_file_struct.session( session_count + 1 ).drug         = drug;
    master_file_struct.session( session_count + 1 ).preprocessed = 0; % Has not yet been preprocessed.
    
    % Save the updated Master Struct if a new struct isn't defined
    if nargin < 6
        save( 'master_file_struct', 'master_file_struct' )
    end
    

end