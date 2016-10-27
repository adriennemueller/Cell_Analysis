function tmp_bar_comp( no_drug, drug )
    ms = 550;
    
    no_drug_avg = mean((sum(no_drug,2) / ms) * 1000);
    no_drug_ste = std((sum(no_drug,2) / ms) * 1000) / sqrt(size(no_drug, 1)); 
    
    drug_avg    = mean((sum(drug,2) / ms) * 1000);
    drug_ste    = std((sum(drug,2) / ms) * 1000) / sqrt(size(drug, 1)); 
    
    
    figure(); bar( [ drug_avg no_drug_avg] ); hold on; 
    errorbar( [drug_avg no_drug_avg], [drug_ste no_drug_ste], '.k');
    
%     vec = [drug_avg no_drug_avg];
%     %fHand = gcf;
%     aHand = axes('parent', gcf);%, fHand);
%     colors = [0 0 0; 1 0 0];
%     for i = 1:numel(vec)
%         bar(i, vec(i), 'parent', aHand, 'facecolor', colors(i,:));
%     end

end