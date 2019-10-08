% Results of this were empty; not bothering to make it spit out filtered
% master file struct.
function signif_attend_mfs = filter_signif_attend( mfs )

    if nargin < 1
        load( 'master_file_struct', 'master_file_struct' );
        mfs = master_file_struct;
    end
    
    signif_attend_mfs = {};
    
    % Run through mfs, getting the filenames, full file paths, and drug of the different processed files
    fn_ffps = get_ffps( master_file_struct );
    fnames = fn_ffps(1,:); ffpaths = fn_ffps(2,:); % drug = fn_ffps(3,:); currents = fn_ffps(4,:);
    
    adjustment_factor = 40; % Num sessions
    
    % Get processed files itemized in Master File Struct, collecting percent
    % correct for cued and uncued trials
    for i = 1:length(fnames)
        
        disp(fnames(i));
        
        % Testing one file
%         if ~ strcmp( fnames(i), 'PROC_2016.08.26_Garfunkel_10mM_SKF81297_3368-4843_MultiUnit.mat')
%             continue
%         end

        load( ffpaths{i}, 'data_struct' ); % Loads data_struct
        
        TrialErrors = [data_struct.trial_error];
        nocue       = [data_struct.nocue];
        
        tes_cue_corrs  = find( (TrialErrors == 0) & (nocue == 0) );
        tes_cue_FPs = find( (TrialErrors == 4) & (nocue == 0) );
        tes_cue_FNs = find( (TrialErrors == 5) & (nocue == 0) );
        n_cue = sum( [length(tes_cue_corrs), length(tes_cue_FPs), length(tes_cue_FNs)]);

        tes_nocue_corrs  = find( (TrialErrors == 0) & (nocue == 1) );
        tes_nocue_FPs = find( (TrialErrors == 4) & (nocue == 1) );
        tes_nocue_FNs = find( (TrialErrors == 5) & (nocue == 1) );
        n_nocue = sum( [length(tes_nocue_corrs), length(tes_nocue_FPs), length(tes_nocue_FNs)]);
    
        perc_Cue = length(tes_cue_corrs) / n_cue;
        perc_NoCue = length(tes_nocue_corrs) / n_nocue; 
    
    
        % Do a chi square test / fisher exact test on the proportions for each
        % file
        pval = chisqcue(length(tes_cue_corrs), (length(tes_cue_FPs) + length(tes_cue_FNs)), length(tes_nocue_corrs), (length(tes_nocue_FPs) + length(tes_nocue_FNs)) );
    
        adjusted_alpha = 0.05 / adjustment_factor;
        
        % Spit out a reduced master_file_stuct, with only the units/proccessed files that are associated
        % with significantly better behavior in cued vs neutrally cued trials
        if pval < adjusted_alpha
            signif_attend_mfs = [signif_attend_mfs fnames(i)];
        end
        
    end

end

% Run through mfs, getting the filenames, full file paths, and drug of the different processed files   
function ffps = get_ffps( mfs )
    
    ffps = {};
    fnames = {};
    drug = {};
    currents = {};
    
    if strcmp( comp_mac_address, 'iMac'), save_direc = '/Users/Adrienne/Documents/MATLAB/Cell_Analysis/Processed';
    else save_direc = '/Users/adrienne/Documents/MATLAB/Cell_Analysis/Processed';
    end

    for i = 1:length(mfs.session)
        sub_direc = mfs.session(i).sub_direc;
        
        for j = 1:length( mfs.session(i).processed_files )
            file =  [mfs.session(i).processed_files(j)];
            fnames = [fnames file];
            
            ffp = fullfile( save_direc, sub_direc, file );
            ffps = [ffps ffp];
            
            curr_drug = mfs.session(i).drug;
            drug = [drug curr_drug];
            
            current = mfs.session(i).currents(j);
            currents = [currents current];
        end
    end
    
    ffps = [fnames; ffps; drug; currents];
end
