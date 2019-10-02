% Wrapper function for gen_overviw_fig - to loop through all processd files
function gen_overview_figs()

    load( 'master_file_struct', 'master_file_struct' );

    % Run through mfs, getting the filenames, full file paths, and drug of the different processed files
    fn_ffps = get_ffps( master_file_struct );
    fnames = fn_ffps(1,:); ffpaths = fn_ffps(2,:); drug = fn_ffps(3,:); currents = fn_ffps(4,:);
   
    % Get list of all unique paradigms in master_file_struct
    %paradigms = horzcat(master_file_struct.session.paradigms);
    %paradigms = unique([paradigms{:}]); % Contains Unknown, so don't want this til solve what Unknown is. Probably probe trials.
    
    % Plot Overview Figs For Each File
    for i = 1:length(fnames)
        
        disp(fnames(i));
        
        % Testing one file
%         if ~ strcmp( fnames(i), 'PROC_2016.05.25_Garfunkel_SKF81297_10mM_6360-7616_MultiUnit1.mat' )
%             continue
%         end

        load( ffpaths{i}, 'data_struct' ); % Loads data_struct
        
        overview_fig = gen_overview_fig( data_struct, currents{i} );
        
        % Save out fig in appropriate directory %
        current_list = currents{i};
        save_name_mat = strcat( 'tmp_figs/',fnames(i), '_', num2str(current_list(2)), 'nA' );
        save_name = strrep(save_name_mat,'.mat','');
        overview_fig.PaperPositionMode = 'auto';
        overview_fig_pos = overview_fig.PaperPosition;
        overview_fig.PaperSize = [overview_fig_pos(3) overview_fig_pos(4)];
        print( overview_fig, strcat(save_name{1}, '.pdf'), '-dpdf'); % .png and .fig also posisble. % May Want to add paradigm to this eventually
        
        close all;
    end

end

% Run through mfs, getting the filenames, full file paths, and drug of the different processed files   
function ffps = get_ffps( mfs )
    
    ffps = {};
    fnames = {};
    drug = {};
    currents = {};
    
    if strcmp( comp_mac_address, 'iMac'), save_direc = '/Users/Adrienne/Documents/MATLAB/Cell_Analysis/Processed';
    else save_direc = '/Users/adrienne/Documents/MATLAB/Cell_Analysis/Processed';
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
