% drug_struct = [sch20, sch50, skf20, skf50];

function plot_attend_svm( drug_struct )


    figure()
  
    for i = 1:4
        ctrl_vec = [1 3 5 7];
        drug_vec = [2 4 6 8];
        
        ctrl_unscram = [drug_struct.avg_perc_corr];
        ctrl_unscram = ctrl_unscram( ctrl_vec(i) );
        ctrl_scram   = [drug_struct.avg_scrambled_perc_corr];
        ctrl_scram   = ctrl_scram( ctrl_vec(i) );
        ctrl_perc_corr = [ctrl_unscram ctrl_scram];
        
        drug_unscram = [drug_struct.avg_perc_corr];
        drug_unscram = drug_unscram( drug_vec(i) );
        drug_scram = [drug_struct.avg_scrambled_perc_corr];
        drug_scram = drug_scram( drug_vec(i) );
        drug_perc_corr = [drug_unscram drug_scram];

               
        subplot( 1,4,i);    
        hold on;
        bar([ctrl_perc_corr; drug_perc_corr] * 100);
              
        %xlabel( 'Contrast' );
        TickLabel_FontSize = 12;
        ylabel( 'Classifier Performance (% Correct)' );
        set( gca, 'YTick', [40 50 60 70], 'XTick', [1, 2], 'XTickLabel', {'Control', 'Drug'}, 'FontSize', TickLabel_FontSize,  ... 
            'FontWeight', 'Bold' );
        xlim([0 3]); 
        ylim([40 70] );
        box( gca, 'off');
        hold off;
        
        
    end



end