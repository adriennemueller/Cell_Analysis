% Wrapper function for gen_overviw_fig - to loop through all processd files
function gen_overview_figs()

    load( 'master_file_struct', 'master_file_struct' );

    % Run through mfs, getting the filenames, full file paths, and drug of the different processed files
    fn_ffps = get_ffps( master_file_struct );
    fnames = fn_ffps(1,:); ffpaths = fn_ffps(2,:); drug = fn_ffps(3,:); currents = fn_ffps(4,:);
   
    % Get list of all unique paradigms in master_file_struct
    
    paradigms = {}; %%% TODO TODO
    
    % Plot Overview Figs For Each File
    for i = 1:length(fnames)
        
        disp(fnames(i));
        load( ffpaths{i} );

        % Loop through all paradigms in that file and gen over view fig
        file_paradigms = {};
        for j = 1:length(file_paradigms)
        
            paradigm = file_paradigms(j);
            overview_fig = gen_overview_fig( data_struct, paradigm )
    
            % Saveout fig in appropriate directory

%             save_name_mat = strcat('tmp_figs/',fnames(i), '_', num2str(attend_struct(j).current));
%             save_name = strrep(save_name_mat,'.mat','');
%               
%             saveas( att_sden_modded_fig, strcat(save_name{1}, '_Modded_DOn.png') );
%             saveas( att_sden_modded_fig, strcat(save_name{1}, '_Modded_DOn.svg') );
%             savefig( att_sden_modded_fig, strcat(save_name{1}, '_Modded_DOn.fig') );              
    
    
    
        end
    end

end

% Run through mfs, getting the filenames, full file paths, and drug of the different processed files   
function ffps = get_ffps( mfs )
    
    ffps = {};
    fnames = {};
    drug = {};
    currents = {};
    
    if strcmp( comp_mac_address, 'iMac'), save_direc = '/Users/Adrienne/Documents/MATLAB/Cell_Analysis/Processed';
    else save_direc = '/Users/eddi/Documents/MATLAB/Cell_Analysis/Processed';
    end

    for i = 1:length(mfs.session)
        sub_direc = mfs.session(i).sub_direc;
        
        for j = 1:length( mfs.session(i).processed_files )
            file =  [mfs.session(i).processed_files(j)];
            fnames = [fnames file];
            
            ffp = fullfile( save_direc, sub_direc, file );
            ffps = [ffps ffp];
            
            curr_drug = mfs.session(i).drug;
            drug = [drug curr_drug];
            
            current = mfs.session(i).currents(j);
            currents = [currents current];
        end
    end
    
    ffps = [fnames; ffps; drug; currents];
end
