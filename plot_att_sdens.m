function att_sden_fig = plot_att_sdens( attend_struct )

    figure('Visible','off');
    hold on;

    n_plots  = length(attend_struct.sden_summs); % Num Directions
    
    % For y-axis height
    max_sden = max( [[attend_struct.sden_summs.dOaI_avg] [attend_struct.sden_summs.dOaO_avg], ...
                    [attend_struct.sden_summs.dNaI_avg] [attend_struct.sden_summs.dNaO_avg]] );

    for i = 1:n_plots

        % Get the drug off/on, attend in/out millis, avg traces and stes
        dOaI_avg = attend_struct.sden_summs(i).dOaI_avg;
        dOaI_ste = attend_struct.sden_summs(i).dOaI_ste;
        dOaI_mil = attend_struct.sden_summs(i).dOaI_millis;

        dOaO_avg = attend_struct.sden_summs(i).dOaO_avg;
        dOaO_ste = attend_struct.sden_summs(i).dOaO_ste;
        dOaO_mil = attend_struct.sden_summs(i).dOaO_millis;

        dNaI_avg = attend_struct.sden_summs(i).dNaI_avg;
        dNaI_ste = attend_struct.sden_summs(i).dNaI_ste;
        dNaI_mil = attend_struct.sden_summs(i).dNaI_millis;

        dNaO_avg = attend_struct.sden_summs(i).dNaO_avg;
        dNaO_ste = attend_struct.sden_summs(i).dNaO_ste;
        dNaO_mil = attend_struct.sden_summs(i).dNaO_millis;
        
        % plot them with the overlay
        subplot(3,3, get_subplotidx(i));
        hold on;
        
        % Delimiting cue window
        line([0 0], [0 max_sden+5], 'LineStyle', ':', 'Color', [0.2 0.2 0.2]);%, 'YLimInclude', 'off', 'XLimInclude', 'off');
        line([500 500], [0 max_sden+5], 'LineStyle', ':', 'Color', [0.2 0.2 0.2]);%, 'YLimInclude', 'off', 'XLimInclude', 'off');
        
        shadedErrorBar(dOaI_mil, dOaI_avg, dOaI_ste, 'k', 1);
        shadedErrorBar(dOaO_mil, dOaO_avg, dOaO_ste, {'--', 'Color', [0.5 0.5 0.5]}, 1);
        shadedErrorBar(dNaI_mil, dNaI_avg, dNaI_ste, 'r', 1);
        shadedErrorBar(dNaO_mil, dNaO_avg, dNaO_ste, {'--', 'Color', [0.5 0 0]}, 1);
        
        ylim([0 max_sden+5] );
        xlim([-500 750]);
        
        set(gca,'Xtick',[-300 0 500],'XTickLabel',{ ['Targets' char(10) 'On'], ['Cue' char(10) 'On'], 'Cue \newline Off' })

        % stick the two d's on the subplot
        dO_dprime = attend_struct.dmat(i,2);
        dN_dprime = attend_struct.dmat(i,3);
        
        text(700,max_sden, num2str(round(dO_dprime,2)));
        text(700,max_sden-5, num2str(round(dN_dprime,2)), 'Color','red');
        
        %h = legend( [], [], 'Drug Off, Attend In', 'Drug Off, Attend Out', 'Drug On, Attend In', 'Drug On, Attend Out');
        
        hold off;
    
    end
    hold off;
    
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