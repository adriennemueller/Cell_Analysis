function stat_struct = process2()

    load('mfs');
    stat_struct = struct;

    % Run through mfs, getting the filenames, full file paths, and drug of the different processed files
    fn_ffps = get_ffps( master_file_struct );
    fnames = fn_ffps(1,:); ffpaths = fn_ffps(2,:); drug = fn_ffps(3,:); currents = fn_ffps(4,:);
   
    
    % Run AttIn_AttOut for that file
    for i = 1:length(fnames)
        
        disp(fnames(i));
        load( ffpaths{i} );

        contrasts_list = {'2017.06.15', '2017.06.18', '2017.06.18', '2017.06.20' };
        contrast_file = find(contains(fnames(i), contrasts_list));
        if contrast_file
           c_attend_struct = attIn_attOut_Contrasts( data_struct, currents{i} );
        end
        attend_struct = attIn_attOut( data_struct, currents{i} ); %% GOOD
                   
        for j = 1:length(attend_struct)
            if isempty( attend_struct(j).sden_summs ), continue, end
            att_sden_modded_fig = plot_att_sdens_Modified(attend_struct(j));%, fnames(i), num2str(attend_struct(j).current));
            
            save_name_mat = strcat('tmp_figs/',fnames(i), '_', num2str(attend_struct(j).current));
            save_name = strrep(save_name_mat,'.mat','');
              
            saveas( att_sden_modded_fig, strcat(save_name{1}, '_Modded_DOn.png') );
            saveas( att_sden_modded_fig, strcat(save_name{1}, '_Modded_DOn.svg') );
            savefig( att_sden_modded_fig, strcat(save_name{1}, '_Modded_DOn.fig') );              

            if contrast_file
                if isempty( c_attend_struct(j).sden_summs ), continue, end
                contrast_fig = plot_att_sdens_Modified_Contrasts(c_attend_struct(j), fnames(i), num2str(c_attend_struct(j).current));
                saveas( contrast_fig, strcat(save_name{1}, '_Modded_Contrasts.png') );
                saveas( contrast_fig, strcat(save_name{1}, '_Modded_Contrasts.svg') );
                savefig( contrast_fig, strcat(save_name{1}, '_Modded_Contrasts.fig') );              
            end
  
            
        end
        

        stat_struct(i).filename = fnames(i);
        stat_struct(i).attend = attend_struct;
        stat_struct(i).drug = drug(i);

        if (isempty(attend_struct)), continue, end %% GOOD
        
        
    end
    
    save( 'stat_struct', 'stat_struct');
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

