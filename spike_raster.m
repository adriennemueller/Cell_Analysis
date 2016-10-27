function rslt = spike_raster( spikes, drug )

    no_drug_mat = [];
    drug_mat = [];
    
    for i = 1:length(spikes)
        aligned_spikes = spikes{i}(800:1350);        
        if drug(i) == -15      
            no_drug_mat = horzcat( no_drug_mat, aligned_spikes);
        else
            drug_mat = horzcat( drug_mat, aligned_spikes);
        end
    end

    no_drug_mat = logical(no_drug_mat'); assignin('base', 'no_drug_mat', no_drug_mat);
    drug_mat = logical(drug_mat'); assignin('base', 'drug_mat', drug_mat);

    
    figure(); 
    subplot(2,1,1); hold on;
    plotSpikeRaster( drug_mat(1:50,:), 'PlotType','vertline' );
    
    subplot(2,1,2); hold on;
    plotSpikeRaster( no_drug_mat(1:50,:), 'PlotType','vertline');
    
    
    %line([0 0], [0 max_sden+5], 'LineStyle', ':', 'Color', [0.2 0.2 0.2]);%,
    
end