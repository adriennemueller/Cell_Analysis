% Loads the chosen processed file
function data_struct = load_processed_file( sub_direc, processed_fname )

    if strcmp( comp_mac_address, 'iMac'), save_direc = '/Users/Adrienne/Documents/MATLAB/Cell_Analysis/Processed';
    else save_direc = '/Users/adrienne/Documents/MATLAB/Cell_Analysis/Processed';
    end
     
    ffp = fullfile( save_direc, sub_direc, processed_fname );
    load( ffp, 'data_struct' );
end