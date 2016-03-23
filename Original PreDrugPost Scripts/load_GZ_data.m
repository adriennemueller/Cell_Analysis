function load_GZ_data()
    stDevinMS = 2000;

    mainpath = '/Users/eddi/Desktop/Lab Meeting 2015.06.12/plx_tmp/';

    celllist = {'Garf_2015.04.22_U001', 'Garf_2015.04.22_MU001', 'Garf_2015.04.22_MU002', ...
                'Garf_2015.04.22_MU003', 'Garf_2015.04.22_MU004_U01', 'Garf_2015.04.22_MU004_U02', ...
                'Garf_2015.04.22_MU005_U01', 'Garf_2015.04.22_MU005_U02', '0421/Garf_2015.04.21_MU001', ...
                '0421/Garf_2015.04.21_MU002', '0421/Garf_2015.04.21_MU003' %, '0421/Garf_2015.04.21_MU004', ...
                %'0421/Garf_2015.04.21_MU005', '0421/Garf_2015.04.21_MU006'
                };
            
    timelist = [ 0 65.3 75.37 100 ; 225 300 319.208 395 ; 1950 1998 2025.9 2075 ; ...
                  6518 6542.3 6590.1 6627; 7009 7033.596 7053.682 7070; 7009 7033.596 7053.682 7070; ...
                  7208.7 7234.02 7250.906 7300; 7208.7 7234.02 7250.906 7300; 1750 1789.636 1811.758 1850; ...
                  2100 2150.860 2178.295 2230; 2250 2290.588 2320.606 2370; 2830 2865.804 2880.272 2915; ...
                  3160 3210.632 3235.816 3280; 3400 3441.081 3469.446 3520
                ];
            
    avgslist = [];
    
    for i = 1:length(celllist)
        
        pregz = load( strcat( mainpath, celllist{i}, '_PreGZ.mat' ) );
        gz = load( strcat( mainpath, celllist{i}, '_GZ.mat' ) );
        postgz = load( strcat( mainpath, celllist{i}, '_PostGZ.mat' ) );
        
        pregz_r = raster(pregz.SPK09(:,1), timelist(i,1), timelist(i,2) );
        gz_r = raster(gz.SPK09(:,1), timelist(i,2), timelist(i,3) );
        postgz_r = raster(postgz.SPK09(:,1), timelist(i,3), timelist(i,4) );
        combined_r = cat(2, pregz_r, gz_r, postgz_r );
        
        pregz_sd = spike_density(pregz_r, stDevinMS);
        gz_sd = spike_density(gz_r, stDevinMS);
        postgz_sd = spike_density(postgz_r, stDevinMS);
        combined_sd = spike_density(combined_r, stDevinMS);
        
        avg_pregz = sum(pregz_r) / (timelist(i,2) - timelist(i,1));
        avg_gz = sum(gz_r) / (timelist(i,3) - timelist(i,2));
        avg_postgz = sum(postgz_r) / (timelist(i,4) - timelist(i,3));
        
        avgslist = [ avgslist ; avg_pregz avg_gz avg_postgz ];
        
        %Plot Spike Overlay
        spikemat = vertcat(pregz.SPK09(:,2:end), gz.SPK09(:,2:end), postgz.SPK09(:,2:end) );
        [f1, a1] = spike_overlay( spikemat );
        
        % Plot Spike Density trace
        [f2, a2] = PreDrugPost_plot(pregz_sd, gz_sd, postgz_sd, combined_sd, celllist{i}); 
        
        
        
        %figure(f1);
                
        hf2 = figure();
        s1 = subplot(2,1,1);
        pos = get(s1, 'Position');
        delete(s1);
        
        hax2=copyobj(a1,hf2);
        set(hax2, 'Position', pos);
        set(gca, 'FontSize', 18);
        
        s2 = subplot(2, 1, 2);
        pos2 = get(s2, 'Position');
        delete(s2);
        
        hax3 = copyobj(a2,hf2);
        set(hax3,'Position',pos2);
        set(gca, 'FontSize', 18);
        
        %copyobj(allchild(a1),ax1);
        %ax2 = subplot(2, 1, 2);
        %copyobj(allchild(a2),ax2);
        close(f1); close(f2);
        


        
    end
    
    plot_avgFR( avgslist );
    
    norm_avgs = normFR( avgslist );
    plot_avgFR( norm_avgs );
    
end

function rslt = normFR( mat )

    [m n] = size(mat);
    
    outmat = ones(m,n);
    
    
    for i = 1:m
       
        outmat(i,2) = mat(i,2) / mat(i,1);
        outmat(i,3) = mat(i,3) / mat(i,1);
        
    end

    rslt = outmat;
    
end












