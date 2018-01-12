function att_sden_fig = plot_att_sdens_Modified_Contrasts( attend_struct, fname, curr_current )

   

    n_directions  = length(attend_struct.sden_summs); % Num Directions
    
    % For y-axis height
    %max_sden = max( [[attend_struct.sden_summs.C25dOaI_avg] [attend_struct.C25sden_summs.dOaO_avg], ...
    %                [attend_struct.sden_summs.C25dNaI_avg] [attend_struct.C25sden_summs.dNaO_avg]] );

    for i = 1:n_directions
        
        figure('Visible','off','rend','painters','pos',[10 10 1000 300]);
        
        hold on;
        
        %%% C10 %%%
        % Get the drug off/on, attend in/out millis, avg traces and stes
        subplot(1,4,1)
        C10dOaI_avg = attend_struct.sden_summs(i).C10dOaI_avg;
        C10dOaI_ste = attend_struct.sden_summs(i).C10dOaI_ste;
        C10dOaI_mil = attend_struct.sden_summs(i).C10dOaI_millis;

        C10dOaO_avg = attend_struct.sden_summs(i).C10dOaO_avg;
        C10dOaO_ste = attend_struct.sden_summs(i).C10dOaO_ste;
        C10dOaO_mil = attend_struct.sden_summs(i).C10dOaO_millis;

        C10dNaI_avg = attend_struct.sden_summs(i).C10dNaI_avg;
        C10dNaI_ste = attend_struct.sden_summs(i).C10dNaI_ste;
        C10dNaI_mil = attend_struct.sden_summs(i).C10dNaI_millis;

        C10dNaO_avg = attend_struct.sden_summs(i).C10dNaO_avg;
        C10dNaO_ste = attend_struct.sden_summs(i).C10dNaO_ste;
        C10dNaO_mil = attend_struct.sden_summs(i).C10dNaO_millis;
        
        
        C10dN_avg = attend_struct.sden_summs(i).C10dN_avg;
        C10dN_ste = attend_struct.sden_summs(i).C10dN_ste;
        C10dN_mil = attend_struct.sden_summs(i).C10dN_millis;

        C10dO_avg = attend_struct.sden_summs(i).C10dO_avg;
        C10dO_ste = attend_struct.sden_summs(i).C10dO_ste;
        C10dO_mil = attend_struct.sden_summs(i).C10dO_millis;
        

        % plot them with the overlay
        %subplot(3,3, get_subplotidx(i));
            hold on;
        
        % Delimiting cue window
        %line([500 500], [0 max_sden+5], 'LineStyle', ':', 'Color', [0.2 0.2 0.2]);%, 'YLimInclude', 'off', 'XLimInclude', 'off');
        %line([1000 1000], [0 max_sden+5], 'LineStyle', ':', 'Color', [0.2 0.2 0.2]);%, 'YLimInclude', 'off', 'XLimInclude', 'off');
        
        % As long as there were enough trials to get an actual error bar,
        % plot an error bar
        if length(C10dOaI_ste)>1 && length(C10dOaO_ste)>1  && length(C10dNaI_ste)>1  && length(C10dNaO_ste)>1 
            plot([1:500], C10dO_avg(1:500), 'k', 'LineWidth',2);
            plot([1:500], C10dN_avg(1:500), 'r', 'LineWidth',2);
            x = [501:1000, 1000:-1:501];        % x, forwards and backwards
            C10attOut = C10dOaO_avg(501:1000);
            C10attIn = C10dOaI_avg(501:1000);
            yy = [C10attOut, fliplr(max(C10attOut, C10attIn))]; % Draw where attIn > attOut
            fill(x,yy,'k', 'FaceAlpha', 0.3, 'LineStyle', 'none');
            yy = [C10attOut, fliplr(min(C10attOut, C10attIn))]; % Draw where attIn < attOut
            fill(x,yy,'k', 'FaceAlpha', 0.1, 'LineStyle', 'none');

            plot([501:1000], C10attOut, 'k:');
            plot([501:1000], C10attIn, 'k');

            C10attNOut = C10dNaO_avg(501:1000);
            C10attNIn = C10dNaI_avg(501:1000);
            yy = [C10attNOut, fliplr(max(C10attNOut, C10attNIn))]; % Draw where attIn > attOut
            fill(x,yy,'r', 'FaceAlpha', 0.3, 'LineStyle', 'none');
            yy = [C10attNOut, fliplr(min(C10attNOut, C10attNIn))]; % Draw where attIn < attOut
            fill(x,yy,'r', 'FaceAlpha', 0.1, 'LineStyle', 'none');
            
            plot([501:1000], C10attNOut, 'r:');
            plot([501:1000], C10attNIn, 'r');
            
        end
        
       

        
%        ylim([0 max_sden+5] );
        %xlim([-500 750]);
        
        % Should maybe use text to make labels centered.
        set(gca,'Xtick',[200 500 1000],'XTickLabel',{ ['Targets \newline On'], ['Cue \newline On'], ['Targ \newline Flip'] });
        

        
        
        %%% C15 %%%
                % Get the drug off/on, attend in/out millis, avg traces and stes
        subplot(1,4,2)
        C15dOaI_avg = attend_struct.sden_summs(i).C15dOaI_avg;
        C15dOaI_ste = attend_struct.sden_summs(i).C15dOaI_ste;
        C15dOaI_mil = attend_struct.sden_summs(i).C15dOaI_millis;

        C15dOaO_avg = attend_struct.sden_summs(i).C15dOaO_avg;
        C15dOaO_ste = attend_struct.sden_summs(i).C15dOaO_ste;
        C15dOaO_mil = attend_struct.sden_summs(i).C15dOaO_millis;

        C15dNaI_avg = attend_struct.sden_summs(i).C15dNaI_avg;
        C15dNaI_ste = attend_struct.sden_summs(i).C15dNaI_ste;
        C15dNaI_mil = attend_struct.sden_summs(i).C15dNaI_millis;

        C15dNaO_avg = attend_struct.sden_summs(i).C15dNaO_avg;
        C15dNaO_ste = attend_struct.sden_summs(i).C15dNaO_ste;
        C15dNaO_mil = attend_struct.sden_summs(i).C15dNaO_millis;
        
        
        C15dN_avg = attend_struct.sden_summs(i).C15dN_avg;
        C15dN_ste = attend_struct.sden_summs(i).C15dN_ste;
        C15dN_mil = attend_struct.sden_summs(i).C15dN_millis;

        C15dO_avg = attend_struct.sden_summs(i).C15dO_avg;
        C15dO_ste = attend_struct.sden_summs(i).C15dO_ste;
        C15dO_mil = attend_struct.sden_summs(i).C15dO_millis;
        

        % plot them with the overlay
        %subplot(3,3, get_subplotidx(i));
        hold on;
        
        % Delimiting cue window
        %line([500 500], [0 max_sden+5], 'LineStyle', ':', 'Color', [0.2 0.2 0.2]);%, 'YLimInclude', 'off', 'XLimInclude', 'off');
        %line([1000 1000], [0 max_sden+5], 'LineStyle', ':', 'Color', [0.2 0.2 0.2]);%, 'YLimInclude', 'off', 'XLimInclude', 'off');
        
        % As long as there were enough trials to get an actual error bar,
        % plot an error bar
        if length(C15dOaI_ste)>1 && length(C15dOaO_ste)>1  && length(C15dNaI_ste)>1  && length(C15dNaO_ste)>1 
            plot([1:500], C15dO_avg(1:500), 'k', 'LineWidth',2);
            plot([1:500], C15dN_avg(1:500), 'r', 'LineWidth',2);
            x = [501:1000, 1000:-1:501];        % x, forwards and backwards
            C15attOut = C15dOaO_avg(501:1000);
            C15attIn = C15dOaI_avg(501:1000);
            yy = [C15attOut, fliplr(max(C15attOut, C15attIn))]; % Draw where attIn > attOut
            fill(x,yy,'k', 'FaceAlpha', 0.3, 'LineStyle', 'none');
            yy = [C15attOut, fliplr(min(C15attOut, C15attIn))]; % Draw where attIn < attOut
            fill(x,yy,'k', 'FaceAlpha', 0.1, 'LineStyle', 'none');

            plot([501:1000], C15attOut, 'k:');
            plot([501:1000], C15attIn, 'k');

            C15attNOut = C15dNaO_avg(501:1000);
            C15attNIn = C15dNaI_avg(501:1000);
            yy = [C15attNOut, fliplr(max(C15attNOut, C15attNIn))]; % Draw where attIn > attOut
            fill(x,yy,'r', 'FaceAlpha', 0.3, 'LineStyle', 'none');
            yy = [C15attNOut, fliplr(min(C15attNOut, C15attNIn))]; % Draw where attIn < attOut
            fill(x,yy,'r', 'FaceAlpha', 0.1, 'LineStyle', 'none');
            
            plot([501:1000], C15attNOut, 'r:');
            plot([501:1000], C15attNIn, 'r');
            
        end
        
%        ylim([0 max_sden+5] );
        %xlim([-500 750]);
        
        % Should maybe use text to make labels centered.
        set(gca,'Xtick',[200 500 1000],'XTickLabel',{ ['Targets \newline On'], ['Cue \newline On'], ['Targ \newline Flip'] });
        
 
        
        %%% C20 %%%
        
                % Get the drug off/on, attend in/out millis, avg traces and stes
        subplot(1,4,3)
        C20dOaI_avg = attend_struct.sden_summs(i).C20dOaI_avg;
        C20dOaI_ste = attend_struct.sden_summs(i).C20dOaI_ste;
        C20dOaI_mil = attend_struct.sden_summs(i).C20dOaI_millis;

        C20dOaO_avg = attend_struct.sden_summs(i).C20dOaO_avg;
        C20dOaO_ste = attend_struct.sden_summs(i).C20dOaO_ste;
        C20dOaO_mil = attend_struct.sden_summs(i).C20dOaO_millis;

        C20dNaI_avg = attend_struct.sden_summs(i).C20dNaI_avg;
        C20dNaI_ste = attend_struct.sden_summs(i).C20dNaI_ste;
        C20dNaI_mil = attend_struct.sden_summs(i).C20dNaI_millis;

        C20dNaO_avg = attend_struct.sden_summs(i).C20dNaO_avg;
        C20dNaO_ste = attend_struct.sden_summs(i).C20dNaO_ste;
        C20dNaO_mil = attend_struct.sden_summs(i).C20dNaO_millis;
        
        
        C20dN_avg = attend_struct.sden_summs(i).C20dN_avg;
        C20dN_ste = attend_struct.sden_summs(i).C20dN_ste;
        C20dN_mil = attend_struct.sden_summs(i).C20dN_millis;

        C20dO_avg = attend_struct.sden_summs(i).C20dO_avg;
        C20dO_ste = attend_struct.sden_summs(i).C20dO_ste;
        C20dO_mil = attend_struct.sden_summs(i).C20dO_millis;
        

        % plot them with the overlay
        %subplot(3,3, get_subplotidx(i));
        hold on;
        
        % Delimiting cue window
        %line([500 500], [0 max_sden+5], 'LineStyle', ':', 'Color', [0.2 0.2 0.2]);%, 'YLimInclude', 'off', 'XLimInclude', 'off');
        %line([1000 1000], [0 max_sden+5], 'LineStyle', ':', 'Color', [0.2 0.2 0.2]);%, 'YLimInclude', 'off', 'XLimInclude', 'off');
        
        % As long as there were enough trials to get an actual error bar,
        % plot an error bar
        if length(C20dOaI_ste)>1 && length(C20dOaO_ste)>1  && length(C20dNaI_ste)>1  && length(C20dNaO_ste)>1 
            plot([1:500], C20dO_avg(1:500), 'k', 'LineWidth',2);
            plot([1:500], C20dN_avg(1:500), 'r', 'LineWidth',2);
            x = [501:1000, 1000:-1:501];        % x, forwards and backwards
            C20attOut = C20dOaO_avg(501:1000);
            C20attIn = C20dOaI_avg(501:1000);
            yy = [C20attOut, fliplr(max(C20attOut, C20attIn))]; % Draw where attIn > attOut
            fill(x,yy,'k', 'FaceAlpha', 0.3, 'LineStyle', 'none');
            yy = [C20attOut, fliplr(min(C20attOut, C20attIn))]; % Draw where attIn < attOut
            fill(x,yy,'k', 'FaceAlpha', 0.1, 'LineStyle', 'none');

            plot([501:1000], C20attOut, 'k:');
            plot([501:1000], C20attIn, 'k');

            C20attNOut = C20dNaO_avg(501:1000);
            C20attNIn = C20dNaI_avg(501:1000);
            yy = [C20attNOut, fliplr(max(C20attNOut, C20attNIn))]; % Draw where attIn > attOut
            fill(x,yy,'r', 'FaceAlpha', 0.3, 'LineStyle', 'none');
            yy = [C20attNOut, fliplr(min(C20attNOut, C20attNIn))]; % Draw where attIn < attOut
            fill(x,yy,'r', 'FaceAlpha', 0.1, 'LineStyle', 'none');
            
            plot([501:1000], C20attNOut, 'r:');
            plot([501:1000], C20attNIn, 'r');
            
        end
        
%        ylim([0 max_sden+5] );
        %xlim([-500 750]);
        
        % Should maybe use text to make labels centered.
        set(gca,'Xtick',[200 500 1000],'XTickLabel',{ ['Targets \newline On'], ['Cue \newline On'], ['Targ \newline Flip'] });
        
        
        
        
        
        %%% C25 %%%
                % Get the drug off/on, attend in/out millis, avg traces and stes
        subplot(1,4,4)
        C25dOaI_avg = attend_struct.sden_summs(i).C25dOaI_avg;
        C25dOaI_ste = attend_struct.sden_summs(i).C25dOaI_ste;
        C25dOaI_mil = attend_struct.sden_summs(i).C25dOaI_millis;

        C25dOaO_avg = attend_struct.sden_summs(i).C25dOaO_avg;
        C25dOaO_ste = attend_struct.sden_summs(i).C25dOaO_ste;
        C25dOaO_mil = attend_struct.sden_summs(i).C25dOaO_millis;

        C25dNaI_avg = attend_struct.sden_summs(i).C25dNaI_avg;
        C25dNaI_ste = attend_struct.sden_summs(i).C25dNaI_ste;
        C25dNaI_mil = attend_struct.sden_summs(i).C25dNaI_millis;

        C25dNaO_avg = attend_struct.sden_summs(i).C25dNaO_avg;
        C25dNaO_ste = attend_struct.sden_summs(i).C25dNaO_ste;
        C25dNaO_mil = attend_struct.sden_summs(i).C25dNaO_millis;
        
        
        C25dN_avg = attend_struct.sden_summs(i).C25dN_avg;
        C25dN_ste = attend_struct.sden_summs(i).C25dN_ste;
        C25dN_mil = attend_struct.sden_summs(i).C25dN_millis;

        C25dO_avg = attend_struct.sden_summs(i).C25dO_avg;
        C25dO_ste = attend_struct.sden_summs(i).C25dO_ste;
        C25dO_mil = attend_struct.sden_summs(i).C25dO_millis;
        

        % plot them with the overlay
        %subplot(3,3, get_subplotidx(i));
        hold on;
        
        % Delimiting cue window
        %line([500 500], [0 max_sden+5], 'LineStyle', ':', 'Color', [0.2 0.2 0.2]);%, 'YLimInclude', 'off', 'XLimInclude', 'off');
        %line([1000 1000], [0 max_sden+5], 'LineStyle', ':', 'Color', [0.2 0.2 0.2]);%, 'YLimInclude', 'off', 'XLimInclude', 'off');
        
        % As long as there were enough trials to get an actual error bar,
        % plot an error bar
        if length(C25dOaI_ste)>1 && length(C25dOaO_ste)>1  && length(C25dNaI_ste)>1  && length(C25dNaO_ste)>1 
            plot([1:500], C25dO_avg(1:500), 'k', 'LineWidth',2);
            plot([1:500], C25dN_avg(1:500), 'r', 'LineWidth',2);
            x = [501:1000, 1000:-1:501];        % x, forwards and backwards
            C25attOut = C25dOaO_avg(501:1000);
            C25attIn = C25dOaI_avg(501:1000);
            yy = [C25attOut, fliplr(max(C25attOut, C25attIn))]; % Draw where attIn > attOut
            fill(x,yy,'k', 'FaceAlpha', 0.3, 'LineStyle', 'none');
            yy = [C25attOut, fliplr(min(C25attOut, C25attIn))]; % Draw where attIn < attOut
            fill(x,yy,'k', 'FaceAlpha', 0.1, 'LineStyle', 'none');

            plot([501:1000], C25attOut, 'k:');
            plot([501:1000], C25attIn, 'k');

            C25attNOut = C25dNaO_avg(501:1000);
            C25attNIn = C25dNaI_avg(501:1000);
            yy = [C25attNOut, fliplr(max(C25attNOut, C25attNIn))]; % Draw where attIn > attOut
            fill(x,yy,'r', 'FaceAlpha', 0.3, 'LineStyle', 'none');
            yy = [C25attNOut, fliplr(min(C25attNOut, C25attNIn))]; % Draw where attIn < attOut
            fill(x,yy,'r', 'FaceAlpha', 0.1, 'LineStyle', 'none');
            
            plot([501:1000], C25attNOut, 'r:');
            plot([501:1000], C25attNIn, 'r');
            
        end
        
%        ylim([0 max_sden+5] );
        %xlim([-500 750]);
        
        % Should maybe use text to make labels centered.
        set(gca,'Xtick',[200 500 1000],'XTickLabel',{ ['Targets \newline On'], ['Cue \newline On'], ['Targ \newline Flip'] });
        
        
        hold off;
    
        
    signif_text = [];
    if sum(attend_struct.anova_mat.p <= 0.05/8) %%% CHECK THIS
        signif_text = 'Significant: ';
        for l = 1:length( attend_struct.anova_mat.p )
            if attend_struct.anova_mat.p(l) <= 0.05/8
                signif_text = strcat(signif_text, {' '}, attend_struct.anova_mat.tbl(l+1,1), {' '}, num2str( attend_struct.anova_mat.p(l) ) );
            end


        end
    end

    if iscell(signif_text)
       annotation('textbox', [0 0.9 1 0.1], 'String', signif_text, 'EdgeColor', 'none', 'HorizontalAlignment', 'center' );
    end    

    att_sden_modded_fig = gcf;
    tightfig(att_sden_modded_fig);
        
    save_name_mat = strcat('tmp_figs/contrast_figs/',fname, '_', num2str(i), '_', curr_current);
    save_name = strrep(save_name_mat,'.mat','');
    saveas( att_sden_modded_fig, strcat(save_name{1}, '_Modded_DOn.png') );
    %savefig( att_sden_modded_fig, strcat(save_name{1}, '_Modded_DOn.fig') );
        
    end
    
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