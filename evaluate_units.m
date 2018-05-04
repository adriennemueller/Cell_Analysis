function rslt = evaluate_units()
    
    % open master file struct
    load( 'tmp_master_file_struct', 'tmp_master_file_struct' );
    
    % if there is no evaluate units field, make one
    if ~ isfield(tmp_master_file_struct.session, 'evaluated')
        tmp_master_file_struct.session(1).evaluated = [];
        [tmp_master_file_struct.session.evaluated] = deal(0);
    end
    
    % find all the sessions that haven't been evaluated
    eval_idxs = find( [tmp_master_file_struct.session.evaluated] == 0 );
    disp(['Unprocessed Sessions = ' num2str( length(eval_idxs) )] );
    
    % go through them and evaluate thems
    for i = 1:length(eval_idxs)
        
        proc_files = tmp_master_file_struct.session(eval_idxs(i)).processed_files;
        
        keep_list = [];
        
        for j = 1:length(proc_files)
            curr_file = proc_files(j);
            curr_fullfile = fullfile( '~/Documents/MATLAB/Cell_Analysis/Processed', tmp_master_file_struct.session(eval_idxs(i)).sub_direc, curr_file);
            
            curr_struct = load(curr_fullfile{1}, 'data_struct');
            curr_struct = curr_struct.data_struct;
            
            currents = tmp_master_file_struct.session(eval_idxs(i)).currents(j);
            
            spike_raster( curr_struct, currents{1} ); drawnow;
            
            %figure();
            %set(0,'CurrentFigure',raster_plot)
            keep = input('Keep this unit? 1/0:','s');
            
            keep_list = [keep_list, str2num(keep)];
        end
        tmp_master_file_struct.session(i).evaluated = 1;
        tmp_master_file_struct.session(i).keep = keep_list;
        save( 'tmp_master_file_struct', 'tmp_master_file_struct' );
    end
end

% function keep = evaluate_unit( unit )
% 
% 
% end