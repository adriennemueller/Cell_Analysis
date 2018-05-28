% This functions generates a master_file struct_from scratch
% Want to update this so it adds new ones without overwriting old master
% file struct, necessarily.

function make_master_file_struct( clean )

    if nargin < 1, clean = 1; end % Clean old mfs by default.

    ionto_super_direc = '/Volumes/Hnoss/Data/Iontophoresis';
    unit_sub_direcs = dir(ionto_super_direc);
    unit_sub_direcs = unit_sub_direcs([unit_sub_direcs.isdir] == 1);
    unit_sub_direcs = unit_sub_direcs( 3:end ); % Eliminate . and .. directories

    % Removes old master_file_struct
    if clean
        disp( 'Cleaning old master_file_struct.');
        delete('master_file_struct.mat');
        master_file_struct = struct;
        master_file_struct.main_direc = '/Users/Adrienne/Documents/MATLAB/Cell_Analysis/';
        save('master_file_struct', 'master_file_struct');
    else
        load( 'master_file_struct.mat', 'master_file_struct' );                
        % If there are folders in master_file_struct that are no longer in
        % the superdirec, remove them from master_file_struct
        mfs_subdirecs = {master_file_struct.session.sub_direc};
        matched_idx_list = ismember( mfs_subdirecs, {unit_sub_direcs.name});
        dud_idxs = find(matched_idx_list == 0);
        master_file_struct.session(dud_idxs) = [];        
    end
        
    
    % Iterate through unit_subdirecs
    for i = 1:length(unit_sub_direcs)
        
        folder = unit_sub_direcs(i);
        search_folder = fullfile( ionto_super_direc, folder.name );
        
        % Identify EventTimes File
        event_file_struct = findfiles( search_folder, '*Event*', 0 );
        if ~ isempty(event_file_struct)
            event_file = event_file_struct.name;    
        else
          disp( strcat('No Event File in: ', {' '}, folder.name ) ); 
          continue
        end
        
        % Identify Unit Files
        % If Sorted2.0 folder exists, search in it
        if isdir( fullfile( search_folder, 'sorted2.0') )
            search_folder = fullfile( search_folder, 'sorted2.0' );
        % Else search in Sorted
        elseif isdir( fullfile( search_folder, 'Sorted') )
            search_folder = fullfile( search_folder, 'Sorted' );
        else
            disp( strcat( 'No Files Sorted in:', {' '}, folder.name ));
            continue
        end
              
        unit_filelist_struct = findfiles( search_folder, '*.mat', 0 );
        unit_filelist = {unit_filelist_struct.name};
        
        % Identify Associated .bhv Files
        folder_date = convert_folder_date(folder.name);
        
        bhv_filelist_struct = findfiles( ionto_super_direc, strcat('*',folder_date,'*'), 0 );
        LP_file_exists = strfind( {bhv_filelist_struct.name}, 'LP' );
        LP_file_indices = find(~cellfun(@isempty, LP_file_exists));
        bhv_filelist_struct = bhv_filelist_struct( LP_file_indices );
        
        if isempty( bhv_filelist_struct )
            disp( strcat('No bhv files in folder:', {' '}, folder.name ) );
            continue;
        end
        bhv_filelist = {bhv_filelist_struct.name};
 
        % Identify Drug
        if ~isempty( strfind( unit_filelist{1}, 'SCH' ) )
            drug = 'SCH23390';
        elseif ~isempty( strfind( unit_filelist{1}, 'SKF' ) )
            drug = 'SKF81297';
        else
            disp( strcat('No Drug Information available in filenames for:', {' '}, folder.name ));
            drug = '';
        end
        
        % Add Session
        master_file_struct = add_session( folder.name, event_file, bhv_filelist, unit_filelist, drug, master_file_struct );
        
    end
    save( 'master_file_struct', 'master_file_struct' );
end

function rslt = convert_folder_date( folder )
    
    year = folder(1:4);
    month = folder(6:7);
    day = folder(9:10);

    rslt = strcat( month, '-', day, '-', year );
end
   
