
function [control_vals, drug_vals, rs_pval] = unity_plot( unity_struct, drug, current, unity_type, position_type )

    control_vals = unity_struct.bestdir_ctrl_vals;
    drug_vals    = unity_struct.bestdir_drug_vals;
    animal = unity_struct.animal;
    signif_visfix = unity_struct.signif_visfix;
    signif_attend = unity_struct.signif_attend;

%     nanidxs = [find(isnan(control_vals) )  find(isnan( drug_vals) ) ];
%      if ~ isempty(nanidxs)
%         rs_control_vals( nanidxs ) = [];
%         rs_drug_vals( nanidxs ) = [];
%      else
%          rs_control_vals = control_vals;
%          rs_drug_vals = drug_vals;
%     end
    
    %rs_pval = ranksum( control_vals, drug_vals );
    rs_pval = signrank( control_vals, drug_vals );
    
    figure('units','normalized','position',[.1 .1 .3 .42]);
    m_color = get_marker_color( drug );
    
    hold on;

    
    %%% Scatter Plot    
    
    % Creat separate lists of Garfunkel and Jose w and w/o specificied
    % significant visfix activity
    garf_idxs = find(strcmp(animal, 'Garfunkel'));
    jose_idxs = find(strcmp(animal, 'Jose')); 
        
    garf_signif_ctrl    = control_vals( intersect( garf_idxs, find(signif_visfix == 1) ) );
    garf_nonsignif_ctrl = control_vals( intersect( garf_idxs, find(signif_visfix ~= 1) ) );
    garf_signif_drug    = drug_vals( intersect(garf_idxs, find(signif_visfix == 1) ) );
    garf_nonsignif_drug = drug_vals( intersect(garf_idxs, find(signif_visfix ~= 1) ) );

    jose_signif_ctrl    = control_vals( intersect(jose_idxs, find(signif_visfix == 1) ) );
    jose_nonsignif_ctrl = control_vals( intersect(jose_idxs, find(signif_visfix ~= 1) ) );
    jose_signif_drug    = drug_vals( intersect(jose_idxs, find(signif_visfix == 1) ) );
    jose_nonsignif_drug = drug_vals( intersect(jose_idxs, find(signif_visfix ~= 1) ) );

    
    
    % Plot all four groups
    plot( garf_signif_ctrl, garf_signif_drug, 'o', 'Color', m_color, 'MarkerFaceColor', m_color );
    plot( garf_nonsignif_ctrl, garf_nonsignif_drug, 'o', 'Color', m_color );

    plot( jose_signif_ctrl, jose_signif_drug, 'd', 'Color', m_color, 'MarkerFaceColor', m_color );
    plot( jose_nonsignif_ctrl, jose_nonsignif_drug, 'd', 'Color', m_color );
    
    % Prettify plot
    ax = gca;
    if strcmp(unity_type, 'd_prime')
        xlabel('Control d'' prime', 'FontSize', 16, 'FontWeight', 'bold'); ylabel( [drug 'd'' prime'], 'FontSize', 16, 'FontWeight', 'bold' );
        xlim([-2 2]); ylim([-2 2]);
        ax.XTick = [-2 -1 0 1 2];
        ax.YTick = ax.XTick;
        line( [-3 3], [-3 3], 'Color', 'black');
        p_x = -1.5; p_y = 1.5;
    elseif strcmp( unity_type, 'mean_fr' )
        xlabel('Control FR', 'FontSize', 16, 'FontWeight', 'bold'); ylabel( [drug 'FR'], 'FontSize', 16, 'FontWeight', 'bold' );
        xlim([0 50]); ylim([0 50]);
        ax.XTick = [0 50];
        ax.YTick = ax.XTick;
        line( [-50 50], [-50 50], 'Color', 'black');
        p_x = 1; p_y = 8;
    elseif strcmp( unity_type, 'mod_idx' )
        xlabel('Control M.I.''', 'FontSize', 16, 'FontWeight', 'bold'); ylabel( [drug 'M.I.'], 'FontSize', 16, 'FontWeight', 'bold' );
        xlim([-1.5 1.5]); ylim([-1.5 1.5]);
        ax.XTick = [-2 -1 0 1 2];
        ax.YTick = ax.XTick;
        line( [-3 3], [-3 3], 'Color', 'black');
        p_x = -1.5; p_y = 1.5;
    end
      
    set(gca,'FontSize',12, 'FontWeight', 'bold');
        
    
    rs_str = get_pval_string(rs_pval);
    ctrl_mean = nanmean( control_vals );
    drug_mean = nanmean( drug_vals );
    
    text(p_x, p_y, ['p ' rs_str], 'FontSize', 16, 'FontWeight', 'bold');
    text(p_x, 1.4, strcat( 'Control Mean = ', num2str(round(ctrl_mean,3)))); 
    text(p_x, 1.3, strcat( 'Drug Mean = ', num2str(round(drug_mean,3)))); 
    text(1.3, -1.7, ['N = ' num2str(length(control_vals))], 'FontSize', 16, 'FontWeight', 'bold');
    hold off;
    title(strcat( drug,{' '} , num2str(current), 'nA' ), 'FontSize', 18, 'Color', m_color);
    
    % Display whether attentional modulation is significant in
    % subpopulation
    [h, garf_p_sig_ctrl] = ttest( garf_signif_ctrl, 0 );
    [h, garf_p_all_ctrl] = ttest( [garf_signif_ctrl, garf_nonsignif_ctrl], 0 );
    [h, garf_p_sig_drug] = ttest( garf_signif_drug, 0 );
    [h, garf_p_all_drug] = ttest( [garf_signif_drug, garf_nonsignif_drug], 0 );

    [h, jose_p_sig_ctrl] = ttest( jose_signif_ctrl, 0 );
    [h, jose_p_all_ctrl] = ttest( [jose_signif_ctrl, jose_nonsignif_ctrl], 0 );
    [h, jose_p_sig_drug] = ttest( jose_signif_drug, 0 );
    [h, jose_p_all_drug] = ttest( [jose_signif_drug, jose_nonsignif_drug], 0 );

    [h, both_p_all_ctrl] = ttest( [garf_signif_ctrl, garf_nonsignif_ctrl, jose_signif_ctrl, jose_nonsignif_ctrl], 0 );
    [h, both_p_all_drug] = ttest( [garf_signif_drug, garf_nonsignif_drug, jose_signif_drug, jose_nonsignif_drug], 0 );

    
    disp( ['Garf Control Signif Vis Mod Units Attend Mod v 0: ', num2str(round(garf_p_sig_ctrl,3))]);
    disp( ['Garf Control All Units Attend Mod v 0: ', num2str(round(garf_p_all_ctrl,3))]);
    disp( ['Garf Drug Signif Vis Mod Units Attend Mod v 0: ', num2str(round(garf_p_sig_drug,3))]);
    disp( ['Garf Drug All Units Attend Mod v 0: ', num2str(round(garf_p_all_drug,3))]);
    g_rs_pval = signrank( garf_signif_ctrl, garf_signif_drug );
    disp( ['Garf Signif Vis Mod Only Sign-Rank p: ', num2str(round(g_rs_pval,3))]);
   
    disp( ['Jose Control Signif Vis Mod Units Attend Mod v 0: ', num2str(round(jose_p_sig_ctrl,3))]);
    disp( ['Jose Control All Units Attend Mod v 0: ', num2str(round(jose_p_all_ctrl,3))]);
    disp( ['Jose Drug Signif Vis Mod Units Attend Mod v 0: ', num2str(round(jose_p_sig_drug,3))]);
    disp( ['Jose Drug All Units Attend Mod v 0: ', num2str(round(jose_p_all_drug,3))]);
    j_rs_pval = signrank( jose_signif_ctrl, jose_signif_drug );
    disp( ['Jose Signif Vis Mod Only Sign-Rank p: ', num2str(round(j_rs_pval,3))]);

    
    disp( ['Both Control All Vis Mod Units Attend Mod v 0: ', num2str(round(both_p_all_ctrl,3))]);
    disp( ['Both Drug All Vis Mod Units Attend Mod v 0: ', num2str(round(both_p_all_drug,3))]);
    b_rs_pval = signrank( [garf_signif_ctrl, jose_signif_ctrl], [garf_signif_drug, jose_signif_drug] );
    disp( ['Both Signif Vis Mod Only Sign-Rank p: ', num2str(round(b_rs_pval,3))]);

    
end

function m_color = get_marker_color( drug )

    if strcmp( drug, 'SCH23390' )
        m_color = 'red';
    elseif strcmp( drug, 'SKF81297' )
        m_color = 'blue';
    end

end


function pval_str = get_pval_string( pval )
    if isempty(pval)
        pval_str = '= NaN';
    elseif pval >= 0.01
        pval_str = [ '= ' num2str(round(pval, 3))];
    elseif (0.001 < pval) && (pval < 0.01)
        pval_str = ['< 0.01'];
    elseif pval < 0.000001
        pval_str = ['< 1*10-6' ];
    elseif pval < 0.001
        pval_str = ['< 0.001'];
    else
        pval_str = '= NaN';
    end
    % Make for smaller
end