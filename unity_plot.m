
function unity_plot( control_vals, drug_vals, drug, current, unity_type )

    rs_pval = ranksum( control_vals, drug_vals );
    
    
    figure('units','normalized','position',[.1 .1 .3 .42]);
    m_color = get_marker_color( drug );
    
    hold on;
    
    %%% Scatter Plot    
    scatterhist( control_vals, drug_vals, 'Color', m_color );
    
    ax = gca;
    if strcmp(unity_type, 'd_prime')
        xlabel('Control d''', 'FontSize', 16, 'FontWeight', 'bold'); ylabel( [drug 'd'''], 'FontSize', 16, 'FontWeight', 'bold' );
        xlim([-2 2]);
        ylim([-2 2]);
        ax.XTick = [-2 -1 0 1 2];
        ax.YTick = ax.XTick;
        set(gca,'FontSize',12, 'FontWeight', 'bold');
        line( [-3 3], [-3 3], 'Color', 'black');
    elseif strcmp( unity_type, 'mean_fr' )
        xlabel('Control FR', 'FontSize', 16, 'FontWeight', 'bold'); ylabel( [drug 'FR'], 'FontSize', 16, 'FontWeight', 'bold' );
        xlim([0 15]);
        ylim([0 15]);
        ax.XTick = [0 10 20 30 40 50];
        ax.YTick = ax.XTick;
        set(gca,'FontSize',12, 'FontWeight', 'bold');
        line( [-50 50], [-50 50], 'Color', 'black');
    end
      
    rs_str = get_pval_string(rs_pval);
    text(-1.5, 1.5, ['p ' rs_str], 'FontSize', 16, 'FontWeight', 'bold');
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
        pval_str = [ '= ' num2str(round(pval, 2))];
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