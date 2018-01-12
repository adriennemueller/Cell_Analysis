function att_sden_fig = plot_att_sdens_Modified( attend_struct )

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
        
        
        dN_avg = attend_struct.sden_summs(i).dN_avg;
        dN_ste = attend_struct.sden_summs(i).dN_ste;
        dN_mil = attend_struct.sden_summs(i).dN_millis;

        dO_avg = attend_struct.sden_summs(i).dO_avg;
        dO_ste = attend_struct.sden_summs(i).dO_ste;
        dO_mil = attend_struct.sden_summs(i).dO_millis;
        

        % plot them with the overlay
        subplot(3,3, get_subplotidx(i));
        hold on;
        
        % Delimiting cue window
        line([500 500], [0 max_sden+5], 'LineStyle', ':', 'Color', [0.2 0.2 0.2]);%, 'YLimInclude', 'off', 'XLimInclude', 'off');
        line([1000 1000], [0 max_sden+5], 'LineStyle', ':', 'Color', [0.2 0.2 0.2]);%, 'YLimInclude', 'off', 'XLimInclude', 'off');
        
        % As long as there were enough trials to get an actual error bar,
        % plot an error bar
        if length(dOaI_ste)>1 && length(dOaO_ste)>1  && length(dNaI_ste)>1  && length(dNaO_ste)>1 
            
            % Plot Visual Part - Since 0 = cue onset; 500 = cue onset, 1000
            % = cue offset
%            shadedErrorBar([1:500], dO_avg(1:500), dO_ste(1:500), 'k', 1);
%            shadedErrorBar([1:500], dN_avg(1:500), dN_ste(1:500), 'r', 1);
            
            plot([1:500], dO_avg(1:500), 'k', 'LineWidth',2);
            plot([1:500], dN_avg(1:500), 'r', 'LineWidth',2);
            
            % Plot Cue Part
%             shadedErrorBar([500:1000], dOaI_avg(500:1000), dOaI_ste(500:1000), 'k', 1);
%             shadedErrorBar([500:1000], dOaO_avg(500:1000), dOaO_ste(500:1000), {'--', 'Color', [0.5 0.5 0.5]}, 1);
%             shadedErrorBar([500:1000], dNaI_avg(500:1000), dNaI_ste(500:1000), 'r', 1);
%             shadedErrorBar([500:1000], dNaO_avg(500:1000), dNaO_ste(500:1000), {'--', 'Color', [0.5 0 0]}, 1);

            x = [501:1000, 1000:-1:501];        % x, forwards and backwards
            
            attOut = dOaO_avg(501:1000);
            attIn = dOaI_avg(501:1000);
            yy = [attOut, fliplr(max(attOut, attIn))]; % Draw where attIn > attOut
            fill(x,yy,'k', 'FaceAlpha', 0.3, 'LineStyle', 'none');
            yy = [attOut, fliplr(min(attOut, attIn))]; % Draw where attIn < attOut
            fill(x,yy,'k', 'FaceAlpha', 0.1, 'LineStyle', 'none');

            plot([501:1000], attOut, 'k:');
            plot([501:1000], attIn, 'k');

            attNOut = dNaO_avg(501:1000);
            attNIn = dNaI_avg(501:1000);
            yy = [attNOut, fliplr(max(attNOut, attNIn))]; % Draw where attIn > attOut
            fill(x,yy,'r', 'FaceAlpha', 0.3, 'LineStyle', 'none');
            yy = [attNOut, fliplr(min(attNOut, attNIn))]; % Draw where attIn < attOut
            fill(x,yy,'r', 'FaceAlpha', 0.1, 'LineStyle', 'none');
            
            plot([501:1000], attNOut, 'r:');
            plot([501:1000], attNIn, 'r');
            
            % lines
%             plot(501:1000, attOut, 'k');
%             plot(501:1000, attIn, 'Color', [0.2 0.2 0.2]);
            
        end
        
        ylim([0 max_sden+5] );
        %xlim([-500 750]);
        
        % Should maybe use text to make labels centered.
        set(gca,'Xtick',[200 500 1000],'XTickLabel',{ ['Targets \newline On'], ['Cue \newline On'], ['Targ \newline Flip'] });

        % stick the two d's on the subplot
        dO_dprime = attend_struct.dmat(i,2);
        dN_dprime = attend_struct.dmat(i,3);
        
        text(600,max_sden, num2str(round(dO_dprime,2)));
        text(600,max_sden-5, num2str(round(dN_dprime,2)), 'Color','red');
        
        % Denote whether cell is visually responsive for this direction
        dO_ranksign = attend_struct.vis_pval(i,2);
        dN_ranksign = attend_struct.vis_pval(i,3);
        title( ['P-Val: DrugOff VR: ' num2str(round(dO_ranksign,2)) ' / DrugOn VR: ' num2str(round(dN_ranksign,2))], 'FontSize', 6 );
        
        %h = legend( [], [], 'Drug Off, Attend In', 'Drug Off, Attend Out', 'Drug On, Attend In', 'Drug On, Attend Out');
        
        hold off;
    end
    
    
     % Plot d's across space
     subplot(3,3, 5);
     hold on;
     plot( attend_struct.dmat(:,2), '-ok' );
     plot( attend_struct.dmat(:,3), '-or' );
     %xlabel( 'Direction' );
     ylabel( 'D''');
     set(get(gca,'YLabel'),'Rotation',0);   
     xlim([0.5 8.5]);
     yl = ylim;
     ylim( [yl(1) (yl(2)+0.2)] );
     set(gca,'Xtick',[1 3 5 7],'XTickLabel',{'Right', 'Up', 'Left', 'Down'});

     % Should put this into a separate function, Also account for
     % multiple comparisons
     sig_val = 0.05;
     drugOff_sig_idxs = find( attend_struct.dmat(:,4) <= sig_val );
     drugOn_sig_idxs  = find( attend_struct.dmat(:,5) <= sig_val );

     text(drugOff_sig_idxs, attend_struct.dmat(drugOff_sig_idxs, 2) + 0.2, '*', 'FontSize', 12);
     text(drugOn_sig_idxs,  attend_struct.dmat(drugOn_sig_idxs,3) + 0.2, '*', 'Color','red', 'FontSize', 12);
        
    
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