function att_sden_fig = plot_att_subbed_sdens( attend_struct )

    figure('Visible','off');
    hold on;

    n_plots  = length(attend_struct.sden_summs); % Num Directions
    
    % For y-axis height
    %max_sden = max( [[attend_struct.sden_summs.dOaI_avg] [attend_struct.sden_summs.dOaO_avg], ...
    %                [attend_struct.sden_summs.dNaI_avg] [attend_struct.sden_summs.dNaO_avg]] );

    for i = 1:n_plots

        % Get the drug off/on, attend in/out millis, avg traces and stes
        dOaI_avg = attend_struct.sden_summs(i).dOaI_avg;
        dOaI_mil = attend_struct.sden_summs(i).dOaI_millis;

        dOaO_avg = attend_struct.sden_summs(i).dOaO_avg;
        dOaO_mil = attend_struct.sden_summs(i).dOaO_millis;

        dNaI_avg = attend_struct.sden_summs(i).dNaI_avg;
        dNaI_mil = attend_struct.sden_summs(i).dNaI_millis;

        dNaO_avg = attend_struct.sden_summs(i).dNaO_avg;
        dNaO_mil = attend_struct.sden_summs(i).dNaO_millis;
        
        
        dOff_subbed = dOaI_avg - dOaO_avg;
        dOn_subbed  = dNaI_avg - dNaO_avg;
        
        ymax = max( max(dOff_subbed), max(dOn_subbed));
        ymin = min( min(dOff_subbed), min(dOn_subbed));
        
        % plot them with the overlay
        subplot(3,3, get_subplotidx(i));
        hold on;
        
        % Delimiting cue window
        line([0 0], [-30 30], 'LineStyle', ':', 'Color', [0.2 0.2 0.2]);%, 'YLimInclude', 'off', 'XLimInclude', 'off');
        line([500 500], [-30 30], 'LineStyle', ':', 'Color', [0.2 0.2 0.2]);%, 'YLimInclude', 'off', 'XLimInclude', 'off');
        

        plot( dOaI_mil, dOff_subbed, 'k');
        plot( dOaI_mil, dOn_subbed, 'r');
        
        %ylim([0 max_sden+5] );
        xlim([-500 750]);
        
        % Should maybe use text to make labels centered.
        set(gca,'Xtick',[-300 0 500],'XTickLabel',{ ['Targets \newline On'], ['Cue \newline On'], ['Cue \newline Off'] });

        % stick the two d's on the subplot
        dO_dprime = attend_struct.dmat(i,2);
        dN_dprime = attend_struct.dmat(i,3);
        
        %text(600,max_sden, num2str(round(dO_dprime,2)));
        %text(600,max_sden-5, num2str(round(dN_dprime,2)), 'Color','red');
        
        % Denote whether cell is visually responsive for this direction
        dO_ranksign = attend_struct.vis_pval(i,2);
        dN_ranksign = attend_struct.vis_pval(i,3);
        title( ['P-Val: DrugOff VR: ' num2str(round(dO_ranksign,2)) ' / DrugOn VR: ' num2str(round(dN_ranksign,2))], 'FontSize', 6 );
        
        %h = legend( [], [], 'Drug Off, Attend In', 'Drug Off, Attend Out', 'Drug On, Attend In', 'Drug On, Attend Out');
        
        hold off;
    end
    
   
    att_sden_fig = gcf;
    
end

% There is a better way
function rslt = get_subplotidx( group )
    if group == 1
        rslt = 6;
    elseif group == 2
        rslt = 3;
    elseif group == 3
        rslt = 2;
    elseif group == 4
        rslt = 1;
    elseif group == 5
        rslt = 4;
    elseif group == 6
        rslt = 7;
    elseif group == 7
        rslt = 8;
    elseif group == 8
        rslt = 9;
    end
end