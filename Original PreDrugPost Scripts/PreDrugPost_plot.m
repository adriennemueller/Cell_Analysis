function [fprop, aprop] = PreDrugPost_plot( pre, drug, post, combined, unitName )

    sden = cat(2, pre, drug, post );

    xs = 1:length(sden);
    xs = xs ./ 1000;
    
    
    pstart = xs(length(pre));
    pend = xs(length(pre) + length(drug));
    pheight = 100;
    
    fprop = figure();
    
    p = patch( [pstart pstart pend pend], [0 pheight pheight 0], 'k' );
    p.FaceAlpha = 0.1;
    p.EdgeColor = 'none';

    hold on;
    

    %plot(xs, sden);
    plot(xs, combined);
    
    ylim( [0 (1.1 *max(sden))] );
    ylabel( 'Spike Density (sp/s)' );
    xlabel( 'Time (s)')
    title( unitName );

    hold off;
    
    aprop = gca;
    
end