function wm_sden_fig = plot_wm_sdens( wm_struct )

    %figure('Visible','off');
    figure();
    hold on;

    n_plots  = length(wm_struct.sden_summs); % Num Directions
    
    % For y-axis height
    max_sden = max( [[wm_struct.sden_summs.wm_avg] [wm_struct.sden_summs.wm_avg], ...
                    [wm_struct.sden_summs.wm_avg] [wm_struct.sden_summs.wm_avg]] );

    for i = 1:n_plots

        % Get the drug off/on, attend in/out millis, avg traces and stes
        wm_avg = wm_struct.sden_summs(i).wm_avg;
        wm_ste = wm_struct.sden_summs(i).wm_ste;
        wm_mil = wm_struct.sden_summs(i).wm_millis;

        % plot them with the overlay
        subplot(3,3, get_subplotidx(i));
        hold on;
        
        % Delimiting cue window
        line([0 0], [0 max_sden+5], 'LineStyle', ':', 'Color', [0.2 0.2 0.2]);%, 'YLimInclude', 'off', 'XLimInclude', 'off');
        line([500 500], [0 max_sden+5], 'LineStyle', ':', 'Color', [0.2 0.2 0.2]);%, 'YLimInclude', 'off', 'XLimInclude', 'off');
        
        % As long as there were enough trials to get an actual error bar,
        % plot an error bar
        if length(wm_ste)>1 
            shadedErrorBar(wm_mil, wm_avg, wm_ste, 'k', 1);

        end
        
        ylim([0 max_sden+5] );
        xlim([-1000 200]);
        
        % Should maybe use text to make labels centered.
        set(gca,'Xtick',[-300 0 500],'XTickLabel',{ ['Targets \newline On'], ['Cue \newline On'], ['Cue \newline Off'] });

        % stick the two d's on the subplot
%        dO_dprime = attend_struct.dmat(i,2);
%        dN_dprime = attend_struct.dmat(i,3);
        
%        text(600,max_sden, num2str(round(dO_dprime,2)));
%        text(600,max_sden-5, num2str(round(dN_dprime,2)), 'Color','red');
        
        % Denote whether cell is visually responsive for this direction
 %       dO_ranksign = attend_struct.vis_pval(i,2);
 %       dN_ranksign = attend_struct.vis_pval(i,3);
 %       title( ['P-Val: DrugOff VR: ' num2str(round(dO_ranksign,2)) ' / DrugOn VR: ' num2str(round(dN_ranksign,2))], 'FontSize', 6 );
        
        %h = legend( [], [], 'Drug Off, Attend In', 'Drug Off, Attend Out', 'Drug On, Attend In', 'Drug On, Attend Out');
        
        hold off;
    end
    
%     
%      % Plot d's across space
%      subplot(3,3, 5);
%      hold on;
%      plot( attend_struct.dmat(:,2), '-ok' );
%      plot( attend_struct.dmat(:,3), '-or' );
%      %xlabel( 'Direction' );
%      ylabel( 'D''');
%      set(get(gca,'YLabel'),'Rotation',0);   
%      xlim([0.5 8.5]);
%      yl = ylim;
%      ylim( [yl(1) (yl(2)+0.2)] );
%      set(gca,'Xtick',[1 3 5 7],'XTickLabel',{'Right', 'Up', 'Left', 'Down'});
% 
%      % Should put this into a separate function, Also account for
%      % multiple comparisons
%      sig_val = 0.05;
%      drugOff_sig_idxs = find( attend_struct.dmat(:,4) <= sig_val );
%      drugOn_sig_idxs  = find( attend_struct.dmat(:,5) <= sig_val );
% 
%      text(drugOff_sig_idxs, attend_struct.dmat(drugOff_sig_idxs, 2) + 0.2, '*', 'FontSize', 12);
%      text(drugOn_sig_idxs,  attend_struct.dmat(drugOn_sig_idxs,3) + 0.2, '*', 'Color','red', 'FontSize', 12);
%         
%     
%     hold off;
%     
    wm_sden_fig = gcf;
    
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