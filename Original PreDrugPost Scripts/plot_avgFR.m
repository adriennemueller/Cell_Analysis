function plot_avgFR( avgslist )

    [m n] = size(avgslist);

    pregz = avgslist( :, 1 );
    gz = avgslist( :, 2 );
    postgz = avgslist( :, 3 );

    avg_pregz = mean(pregz);
    avg_gz = mean(gz);
    avg_postgz = mean(postgz);
    
    figure();
    hold on;
 
    bg = bar( [avg_pregz, avg_gz, avg_postgz] );
    bg.FaceColor = [0.9 0.9 0.9];
     
    for i = 1:m
        
        plot(avgslist(i,:), 'ko-');

    end
    
    ylabel( 'Firing Rate (sp/s)' );
    set(gca,'XTick',[1 2 3]);
    set(gca,'XTickLabel', {'Pre Drug', 'Gabazine', 'Post Drug'})
    set(gca, 'FontSize', 18);
    
    hold off;

end