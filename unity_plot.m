
function [control_vals, drug_vals, rs_pval] = unity_plot( unity_struct, drug, current, unity_type, position_type )

    control_vals = unity_struct.bestdir_ctrl_vals;
    drug_vals    = unity_struct.bestdir_drug_vals;
    animal = unity_struct.animal;
    signif_visfix = unity_struct.signif_visfix;
    signif_attent = unity_struct.signif_attend;

    nanidxs = [find(isnan(control_vals) )  find(isnan( drug_vals) ) ];
     if ~ isempty(nanidxs)
        control_vals( nanidxs ) = [];
        drug_vals( nanidxs ) = [];
    end
    
    rs_pval = ranksum( control_vals, drug_vals );
    
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

    jose_signif_drug    = control_vals( intersect(jose_idxs, find(signif_visfix == 1) ) );
    jose_nonsignif_drug = control_vals( intersect(jose_idxs, find(signif_visfix ~= 1) ) );
    jose_signif_drug    = drug_vals( intersect(jose_idxs, find(signif_visfix == 1) ) );
    jose_nonsignif_drug = drug_vals( intersect(jose_idxs, find(signif_visfix ~= 1) ) );

    % Plot all four groups
    plot( garf_signif_ctrl, garf_signif_drug, 'o', 'Color', m_color );
    plot( garf_nonsignif_ctrl, garf_nonsignif_drug, 'o', 'Color', m_color, 'MarkerFaceColor', m_color );

    plot( jose_signif_ctrl, jose_signif_drug, 's', 'Color', m_color );
    plot( jose_nonsignif_ctrl, jose_nonsignif_drug, 's', 'Color', m_color, 'MarkerFaceColor', m_color );
    
    
    
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
        xlim([-2 2]); ylim([-2 2]);
        ax.XTick = [-2 -1 0 1 2];
        ax.YTick = ax.XTick;
        line( [-3 3], [-3 3], 'Color', 'black');
        p_x = -1.5; p_y = 1.5;
    end
      
    set(gca,'FontSize',12, 'FontWeight', 'bold');
        
    
    
    rs_str = get_pval_string(rs_pval);
    %control_vals( ~isfinite( control_vals )) = NaN;
    ctrl_mean = nanmean( control_vals );
    %drug_vals( ~isfinite( drug_vals )) = NaN;
    drug_mean = nanmean( drug_vals );
    
    text(p_x, p_y, ['p ' rs_str], 'FontSize', 16, 'FontWeight', 'bold');
    text(p_x, 1.4, strcat( 'Control Mean = ', num2str(ctrl_mean))); 
    text(p_x, 1.3, strcat( 'Drug Mean = ', num2str(drug_mean))); 
    text(1.3, -1.7, ['N = ' num2str(length(control_vals))], 'FontSize', 16, 'FontWeight', 'bold');
    hold off;
    title(strcat( drug,{' '} , num2str(current), 'nA' ), 'FontSize', 18, 'Color', m_color);
    
    

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