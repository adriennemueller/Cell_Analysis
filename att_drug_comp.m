function rslt = att_drug_comp( stat_struct, current )


    SCH_drugoff_attmods = [];
    SCH_drugon_attmods  = [];
    
    SKF_drugoff_attmods = [];
    SKF_drugon_attmods  = [];
    
    for i = 1:length(stat_struct )
        
        if isfield(stat_struct(i).attend.anova, 'n_spikes')
        
            disp(stat_struct(i).filename)
            
            if length(stat_struct(i).attend) > 1
                for j = 1:length(stat_struct(i).attend)
                    if (stat_struct(i).attend(j).current == selected_current)
                        attmod_subset = get_attmod_subset( stat_struct(i).attend(j).anova.n_spikes, stat_struct(i).attend(j).dmat );
                        curr_Off_attmods = attmod_subset(:,1);
                        curr_On_attmods  = attmod_subset(:,2);
                    else
                        curr_Off_attmods = [];
                        curr_On_attmods = [];
                    end
                end
            
            else
                    if (stat_struct(i).attend.current == selected_current)
                        attmod_subset = get_attmod_subset( stat_struct(i).attend.anova.n_spikes, stat_struct(i).attend.dmat );
                        curr_Off_attmod = attmod_subset(:,1);
                        curr_On_attmod  = attmod_subset(:,2);
                    else
                        curr_Off_attmod = [];
                        curr_On_attmod = [];
                    end
            
            end

            if strcmp(stat_struct(i).drug, 'SCH23390')
                SCH_drugoff_attmods = [SCH_drugoff_attmods curr_Off_attmod'];
                SCH_drugon_attmods  = [SCH_drugon_attmods  curr_On_attmod' ];
            elseif strcmp(stat_struct(i).drug, 'SKF81297')
                SKF_drugoff_attmods = [SKF_drugoff_attmods curr_Off_attmod'];
                SKF_drugon_attmods  = [SKF_drugon_attmods  curr_On_attmod' ];
            else
                disp('Drug not recognized.')
            end
        end

    end
    
    % Get Rank Sum p-value
    if ~ isempty( SCH_drugoff_dprimes)
        SCH_rs = ranksum( SCH_drugoff_attmods, SCH_drugon_att_mods );
    else
        SCH_rs = [];
    end
    
    if ~ isempty( SKF_drugoff_dprimes)
        SKF_rs = ranksum( SKF_drugoff_attmods, SKF_drugon_att_mods );
    else
        SKF_rs = [];
    end
    
    
    rslt = [SCH_rs SKF_rs];
    %figure();






end

% This function gets the subset of d's which are likely to be the 'attend
% in' d' condition and its two closest neighbours. One point is therefore
% being dropped.
function dprime_subset = get_attmod_subset( n_spikes, dmat )


%%% If Taking All Four
%     Offs = vispval(:,4);
%     Ons  = vispval(:,5);
%     
%     idxs = 1:4;
%     
%     OffVals = Offs(idxs);
%     OnVals  = Ons( idxs );

    drugOn_n_spikes = mean( horzcat(n_spikes.drugOff_AttIn_n_spikes, n_spikes.drugOff_AttOut_n_spikes));
    drugOff_n_spikes = mean( horzcat(n_spikes.drugOff_AttIn_n_spikes, n_spikes.drugOff_AttOut_n_spikes));

    Offs = dmat(:,2);
    Ons  = dmat(:,3);
    
    maxOff_idx = find(Offs == max(Offs));
    if maxOff_idx == 8
        idxs = [7 8 1];
    elseif maxOff_idx == 1
        idxs = [8 1 2];
    else
        idxs = [ (maxOff_idx -1) maxOff_idx (maxOff_idx +1) ];
    end
    
    OffVals = vis_Offs(idxs);
    OnVals  = vis_Ons( idxs );

    dprime_subset = [OffVals, OnVals];

end
