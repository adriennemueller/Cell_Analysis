%%% TODO MAKE THIS INTO a function that processes for different things,
%%% instead of for all things?

% Currently when this is run it will 1) generate d' for drug on/ drug off
% attend in / attend out conditions. and 2) generate cue aligned spike
% density functions for the 4 types. (For all 8 directions);

%%% SHOULD THIS STRUCT also contain the structs for all units, so it can
%%% all be loaded at once? Then, at the end of preprocess, the unit-struct
%%% would be saved to the master_stat_struc (master data struct?) as well.
%%% Also, units could have unique identifiers, instead of just filename.
function stat_struct = process()

    % Load master-file-struct and statistics-struct
    load('master_file_struct');
    load('stat_struct');
        
    % Run through mfs, getting the filenames, full file paths, and drug of the different processed files
    fn_ffps = get_ffps( master_file_struct );
    fnames = fn_ffps(1,:); ffpaths = fn_ffps(2,:); drug = fn_ffps(3,:); currents = fn_ffps(4,:);
    
    %ss_filenames = [stat_struct.filename];
    
    % If any of the filenames don't exist yet (OR, better, if they haven't been
    % attIn/attOut compared yet) - load in that unit's file.

    % Run AttIn_AttOut for that file
    for i = 1:length(fnames)
        
        load( ffpaths{i} );
        attend_struct = attIn_attOut( data_struct, currents{i} );
        
        stat_struct(i).filename = fnames(i);
        stat_struct(i).attend = attend_struct;
        stat_struct(i).drug = drug(i);
       
        if (isempty(attend_struct)), continue, end
        
        % Go through each current
        for j = 1:length(attend_struct)
            if (isempty(find(isinf(attend_struct(j).dmat)))) || ...
                (isempty(find(isnan(attend_struct(j).dmat))))
            att_sden_fig = plot_att_sdens( attend_struct(j) );

            save_name_mat = strcat('tmp_figs/',fnames(i), '_', num2str(attend_struct(j).current));
            save_name = strrep(save_name_mat,'.mat','');
            saveas( att_sden_fig, strcat(save_name{1}, '.png') );
            savefig( att_sden_fig, strcat(save_name{1}, '.fig') );
            end
        end

    end

    % Append the result (d' matrix) to a the summary statistic struct and
    % re-save it out
    save('stat_struct', 'stat_struct');

end


function ffps = get_ffps( mfs )
    
    ffps = {};
    fnames = {};
    drug = {};
    currents = {};
    main_direc = mfs.main_direc;

    for i = 1:length(mfs.session)
        sub_direc = mfs.session(i).sub_direc;
        
        for j = 1:length( mfs.session(i).processed_files )
            file =  [mfs.session(i).processed_files(j)];
            fnames = [fnames file];
            
            ffp = strcat( main_direc, sub_direc, file );
            ffps = [ffps ffp];
            
            curr_drug = mfs.session(i).drug;
            drug = [drug curr_drug];
            
            current = mfs.session(i).currents(j);
            currents = [currents current];
        end
        
    end
    
    ffps = [fnames; ffps; drug; currents];
end