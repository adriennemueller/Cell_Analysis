function svm_plot = plot_svm_drug_effect_histos( epoch_struct )

    epoch_struct_fields = fieldnames(epoch_struct);

    svm_plot = figure();
    svm_plot_pos = get(svm_plot, 'Position'); %// gives x left, y bottom, width, height
    svm_plot_width = svm_plot_pos(3) * length(epoch_struct_fields);  svm_plot_height = svm_plot_pos(4);
    set(svm_plot, 'Position', [10 10 svm_plot_width, svm_plot_height]);

    
    for i = 1:length(epoch_struct_fields)
        
        diff_list = []; % Reset diff_list;
        
        curr_epoch = epoch_struct_fields{i};
        parameter_fields = fieldnames(epoch_struct.(curr_epoch));
        for j = 1:length(parameter_fields)
        
            curr_parameter = parameter_fields{j};
            ctrl_vals = epoch_struct.(curr_epoch).(curr_parameter).ctrl_avg_perc_correct;
            drug_vals = epoch_struct.(curr_epoch).(curr_parameter).drug_avg_perc_correct;
            
            diff_list(j,:) = drug_vals - ctrl_vals;
            
        end
            
        % Once have gone through all cells; have a distribution of performances and
        % compare to chance. (With e.g. an anova? t-test) i.e. see if
        % distribution shifted rightward from 0 mean?
        
        % Plot histos of differences in performances of svm drug vs control
        % with p value of whether center different from 0? Check with Basti.
        subplot( 1, length(epoch_struct_fields), i )
        
        for j = 1:size(diff_list,1)
            
            offset = (j-1);%*2;
            plot_vals = diff_list(j,:) + offset;
            [~, ttestp] = ttest(diff_list(j,:));

            histogram( plot_vals ); 
            hold on;

            yl = ylim;
            new_ymax = yl(2)*1.1;
            %ylim([yl(1), new_ymax]);
            text(offset, yl(2), {'p =' round(ttestp,3)}, 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
            
            % Add vertical lines at "zeros";
            line( [offset offset], [0 100], 'Color', 'black');
            ylim([yl(1), new_ymax]);

        end
        
        xticks( 0 : 1 : (length(parameter_fields)-1) );
        xticklabels( strrep(parameter_fields, '_', ' ') );
        curr_axis = get(gca,'XTickLabel');  
        set(gca,'XTickLabel', curr_axis, 'FontSize', 12, 'FontWeight', 'bold');
        
        yticks(0:5:100);
        ylim([yl(1), new_ymax]); % Reset Yaxis height

        
        if (i == 1) % First Subplot
            ylabel('Number of Units');
        end
            
        % Epoch as title
        title( strrep(curr_epoch, '_', ' ') ) ;
        
    
    end
    
    % Super title with drug and current?
    % sgtitile()
    % title(strcat( drug,{' '} , num2str(current), 'nA' ), 'FontSize', 18, 'Color', m_color);
    
    
end