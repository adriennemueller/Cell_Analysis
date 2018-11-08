function make_drug_plots()

    load( 'sfn_drug_data' );
    
    %sch20_ctrl = sch20_ctrl * -1;
    
    sch20_ctrl_mean = mean(sch20_ctrl);
    sch20_ctrl_std  = std(sch20_ctrl) ./ sqrt(length(sch20_ctrl));
    sch20_drug_mean = mean(sch20_drug);
    sch20_drug_std  = std(sch20_drug) ./ sqrt(length(sch20_drug));

    
    sch50_ctrl_mean = mean(sch50_ctrl);
    sch50_ctrl_std  = std(sch50_ctrl) ./ sqrt(length(sch50_ctrl));
    sch50_drug_mean = mean(sch50_drug);
    sch50_drug_std  = std(sch50_drug) ./ sqrt(length(sch50_drug));
    
    
    skf20_ctrl_mean = mean(skf20_ctrl);
    skf20_ctrl_std  = std(skf20_ctrl) ./ sqrt(length(skf20_ctrl));
    skf20_drug_mean = mean(skf20_drug);
    skf20_drug_std  = std(skf20_drug) ./ sqrt(length(skf20_drug));

    sch20_diff  = sch20_drug - sch20_ctrl;
    sch20_diff_mean = mean( sch20_diff );
    sch20_diff_ste = std( sch20_diff ) ./ sqrt(length(sch20_diff));
    
    sch50_diff  = sch50_drug - sch50_ctrl;
    sch50_diff_mean = mean( sch50_diff );
    sch50_diff_ste = std( sch50_diff ) ./ sqrt(length(sch50_diff));
    
    rs_val = ranksum( sch20_diff, sch50_diff );
    
    %%% Bar Plot
    figure();
    gray_col = [0.8 0.8 0.8];
    
    ctrl = [sch20_ctrl_mean, sch50_ctrl_mean]; ctrl_ste = [sch20_ctrl_std, sch50_ctrl_std];
    drug = [sch20_drug_mean, sch50_drug_mean]; drug_ste = [sch20_drug_std, sch50_drug_std];
    
    b = bar( [0.75, 1.25], ctrl, 'k');
    hold on;
    c = bar( [2.75 3.25], drug, 'r');
 
%     err_x1 = b(1).XData + b(1).XOffset;
%     err_x2 = b(2).XData + b(2).XOffset;
     errorbar( [0.75, 1.25], ctrl, ctrl_ste, 'k', 'linestyle', 'none' );
     errorbar([2.75 3.25], drug, drug_ste, 'k', 'linestyle', 'none' );
%     hold off;

    %xpos1 = [b(1).XData(1) + b(1).XOffset, b(1).XData(1) + b(2).XOffset];
    %xpos2 = [b(1).XData(2) + b(1).XOffset, b(1).XData(2) + b(2).XOffset];
    %sigstar( {xpos1, xpos2}, [D1R_SMI32_pval, D1R_NRG_pval]); 
    
    TickLabel_FontSize = 16; Axis_Font_Size = 20;
    ylabel( 'Attentional Modulation (d'')' );
    set( gca, 'YTick', [-0.1 0 0.1 0.2 0.3], 'XTick', [0.75, 1.25, 2.75, 3.25], 'XTickLabel', {'20nA', '50nA', '20nA', '50nA'}, 'FontSize', TickLabel_FontSize,  ... 
           'FontWeight', 'Bold' );

%    text(2, 0.4, strcat( 'p = ', {' '}, num2str(round(rs_val,3)) ), 'FontSize', Axis_Font_Size, 'FontWeight', 'bold');
    
    %xlim([]); 
   ylim([-0.15 0.3] );

    
    hold off;
    
    %%% Diff Line Plot
    figure();
    plot( [sch20_diff_mean, sch50_diff_mean], 'k-' ); 
    hold on;
    errorbar([1,2], [sch20_diff_mean, sch50_diff_mean], [sch20_diff_ste, sch50_diff_ste], 'k')

    TickLabel_FontSize = 16; Axis_Font_Size = 20;
    ylabel( 'Difference in Attentional Modulation' );
    set( gca, 'YTick', [-0.4 -0.2 0 0.2 0.4], 'XTick', [1, 2], 'XTickLabel', {'20nA', '50nA'}, 'FontSize', TickLabel_FontSize,  ... 
            'FontWeight', 'Bold' );

    text(2, 0.4, strcat( 'p = ', {' '}, num2str(round(rs_val,3)) ), 'FontSize', Axis_Font_Size, 'FontWeight', 'bold');
    
    xlim([0.5 2.5]); 
    ylim([-0.5 0.5] );
    xL = get(gca, 'XLim');
    plot(xL, [0 0], '--k')
    box( gca, 'off');
    hold off;
    
end