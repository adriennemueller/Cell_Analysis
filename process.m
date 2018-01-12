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
        
        %tmp = strfind( fnames(i), '2017.06.18' );
        %if ~isempty( tmp{1})  %%% DO ONLY THIS FOLDER 
        
        disp(fnames(i));
        
        load( ffpaths{i} );
        
%       wm_struct = wmIn_wmOut( data_struct, currents{i} );
%        attend_struct = attIn_attOut( data_struct, currents{i} ); %% GOOD
        attend_struct = attIn_attOut_Contrasts( data_struct, currents{i} );
        for j = 1:length(attend_struct)
            if isempty( attend_struct(j).sden_summs ), continue, end
            plot_att_sdens_Modified_Contrasts(attend_struct(j), fnames(i), num2str(attend_struct(j).current));
        end

        
 %%% UNCOMMENT ONCE FROM HERE TO RETURN 
% %         if (isempty(wm_struct))
% %             continue
% %         end
% %  
% %         % IF NO WM TRIALS
% %         [a b] = size(wm_struct(1).sdens(1).wm);
% %         if (a == 1)
% %             continue
% %         end
% %  
%  %%% GOOD %%%
%        stat_struct(i).filename = fnames(i);
%        stat_struct(i).attend = attend_struct;
%        stat_struct(i).drug = drug(i);
%        
%        if (isempty(attend_struct)), continue, end %% GOOD
% 
% %        wm_sden_fig = plot_wm_sdens( wm_struct );
% %         for j = 1:length(wm_struct)
% %             wm_sden_fig = plot_wm_sdens( wm_struct(j) );
% %             save_name_mat = strcat('tmp_figs/',fnames(i), '_wm_', num2str(wm_struct(j).current));
% %             save_name = strrep(save_name_mat,'.mat','');
% % 
% %             saveas( wm_sden_fig, strcat(save_name{1}, '.png') );
% %             savefig( wm_sden_fig, strcat(save_name{1}, '.fig') );
%             
%             
%             
%        %%% THIS IS GOOD ATTEND CODE  - COMMENTED WHILE HACKED WM  
% %        Go through each current
%         for j = 1:length(attend_struct)
%             if (isempty(find(isinf(attend_struct(j).dmat)))) || ...
%                 (isempty(find(isnan(attend_struct(j).dmat))))
%             
%                        
%                att_sden_fig = plot_att_sdens( attend_struct(j) );
%                att_sden_subbed_fig = plot_att_subbed_sdens( attend_struct(j) );
%                att_sden_modded_fig = plot_att_sdens_Modified(attend_struct(j));
%              
%              % Remove NEW from this after check new attend window.
%              save_name_mat = strcat('tmp_figs/',fnames(i), '_', num2str(attend_struct(j).current));
%              save_name = strrep(save_name_mat,'.mat','');
% %              
% %              saveas( att_sden_fig, strcat(save_name{1}, '.png') );
% %              saveas( att_sden_subbed_fig, strcat(save_name{1}, '_Subbed.png') );
%                saveas( att_sden_modded_fig, strcat(save_name{1}, '_Modded_DOn.png') );
% %              
% %              savefig( att_sden_fig, strcat(save_name{1}, '.fig') );
% %              savefig( att_sden_subbed_fig, strcat(save_name{1}, '_Subbed.fig') );
%                savefig( att_sden_modded_fig, strcat(save_name{1}, '_Modded_DOn.fig') );
% %              
% %              writecsv(attend_struct(j).anova_mat.tbl, save_name{1});
%             end
%             
%          end
       
     %    end
     end
% 
%     % Append the result (d' matrix) to a the summary statistic struct and
%     % re-save it out
 %    save('stat_struct', 'stat_struct'); %%% GOOD

end


function writecsv( tbl, savename )
    filename = strcat(savename, '_Anova.csv');
    fid = fopen(filename, 'w') ;
%     fprintf(fid, '%s,', tbl{1,1:end-1});
%     fprintf(fid, '%s\n', tbl{1,end});
%     
%     fprintf(fid,'%s, %f, %f, %f, %f, %f, %f\n',tbl{2:end,:});

%fprintf('%8s %8s %8s %8s\n', 'col1', 'col2', 'col3', 'col4');
for i = 1:size(tbl,1)
  temp = cellfun(@(x) num2str(x,'%8.4g'),tbl(i,:),'UniformOutput',false);
  fprintf(fid, '%s, %s, %s, %s, %s, %s, %s\n',temp{:});
end



    fclose(fid) ;
     %dlmwrite( filename, tbl(2:end,:), '-append') ;
end


function ffps = get_ffps( mfs )
    
    ffps = {};
    fnames = {};
    drug = {};
    currents = {};
    save_direc = '/Users/eddi/Documents/MATLAB/Cell_Analysis/Processed';

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