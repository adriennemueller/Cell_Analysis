function rslt = corr_att_wm( mfs )

    att_list = [];
    wm_list = [];

    for i = 1:length(mfs.session)

        if isempty( mfs.session(i).attend_stats ) || isempty( mfs.session(i).wm_stats )
                continue
        end      
        
        %num_att_stats = length( mfs.session(i).attend_stats );
        %num_wm_stats  = length( mfs.session(i).wm_stats    );
        
        %num_compare = min(num_att_stats, num_wm_stats);
        
        good_at_idxs = find( ~cellfun(@isempty, mfs.session(i).attend_stats ) ); 
        good_wm_idxs = find( ~cellfun(@isempty,mfs.session(i).wm_stats ) ); 
        
        for j = 1:length( good_wm_idxs )
            wm_idx = good_wm_idxs( j );
            [ d, closest_att_idx ] = min( abs( good_at_idxs - wm_idx ) );
            att_idx = good_at_idxs( closest_att_idx );
            
            att_stats = mfs.session(i).attend_stats{att_idx};
            wm_stats  = mfs.session(i).wm_stats{wm_idx};
            
        
%        for j = 1:num_compare

            
%             att_stats = mfs.session(i).attend_stats{j};
%             wm_stats = mfs.session(i).wm_stats{j};
%             
%             if isempty( att_stats ) || isempty( wm_stats )
%                 continue
%             end      

                for k = 1:length(att_stats)
                    full_ctrl_att_dprimes = [att_stats(k).control_dmat.dmat.dprime_val];
                    ctrl_att_dprimes = full_ctrl_att_dprimes(1:4);
                    
                    full_drug_att_dprimes = [att_stats(k).drug_dmat.dmat.dprime_val];
                    drug_att_dprimes = full_drug_att_dprimes(1:4);
                    
                    full_ctrl_wm_dprimes = [wm_stats(k).control_dmat.dmat.dprime_val];
                    ctrl_wm_dprimes = full_ctrl_wm_dprimes(1:4);
                    
                    full_drug_wm_dprimes = [wm_stats(k).drug_dmat.dmat.dprime_val];
                    drug_wm_dprimes = full_drug_wm_dprimes(1:4);
                    
                    att_vals = drug_att_dprimes - ctrl_att_dprimes;
                    wm_vals = drug_wm_dprimes - ctrl_wm_dprimes;
                    
                    att_list = [att_list att_vals];
                    wm_list  = [wm_list  wm_vals];
                end
                
        end            
            
 %       end
    end

    % Remove outliers
    bad_att_idxs = find( abs(att_list) > 2);
    bad_wm_idxs = find( abs(wm_list) > 10);
    bad_wm_idxs_also = find( isnan( wm_list ) );
    
    bad_idxs = unique( [bad_att_idxs bad_wm_idxs bad_wm_idxs_also] );
    att_list = att_list(:,setdiff(1:end,bad_idxs));
    wm_list  = wm_list(:,setdiff(1:end,bad_idxs));
    
    mdl = fitlm( att_list, wm_list );
    R_Ordinary = mdl.Rsquared.Ordinary;
    mdl_intercept = mdl.Coefficients.Estimate(1);
    mdl_slope = mdl.Coefficients.Estimate(2);
    
    mdl_intercept_p = mdl.Coefficients.pValue(1);
    mdl_slope_p = mdl.Coefficients.pValue(2);
    
    [R_corr pval ] = corrcoef( att_list, wm_list );
    R_corr = R_corr(1,2);
    pval = pval( 1,2);
    
    figure(); plot( att_list, wm_list, 'ok', 'MarkerFaceColor', 'k' );
    hold on;
    
    fit_x = (-2:2);
    fit_y = mdl_slope * fit_x + mdl_intercept;
   
    plot( fit_x, fit_y, '-r' );
    %ylim( [0 120] );
    set(gca, 'ytick', [-6:3:8]);
    set(gca, 'xtick', [-2:1:2]);
    %set(gca,'XTickLabel',{'SMI-32','Neurogranin'});
    %set(gca,'XTickLabel',{'Neurogranin','SMI-32'});
    xlabel( 'Attention d''', 'FontSize', 20, 'FontWeight', 'bold' );
    ylabel( 'Working Memory d''', 'FontSize', 20, 'FontWeight', 'bold' );
    set(gca,'FontSize',14, 'FontWeight', 'bold');
     
    
end