
%%% Want function to make summary anova fig - how many cells with attend
%%% effects exhibit drug/attend interactions. How many drug effects etc.
% Merge gen_anova_counts and gen_anova_contrast_counts


% Should combined drug_svm? (jackknife also)
% MAKE SVM MORE GENERIC FOR other label types (attin/attout, e.g.)
% Make attend_svm not a thing (a hacked copy-paste of drug_svm )
% Make a better version of make_drug_plots for output of SVM analysis.

% Make a better version of plot_contrast_dprimes; which doesn't require
% manual rotation / does something better than the rotation.

% Perform SVM on drug-no drug for different epochs

% Consider unrewarded trials?
Fix sden so it handles spacer gaps.

%%% THINK %%%

% Plot SVM performance as increase number of neurons included

% Plot SVM performance vs scrambled SVM Performance (unity plot? or two
% histograms)

WM DOES NOT HAVE 8 UNIQUE VALUES FOR EPOCHS. FIX.
Make anova also hand attend direction as a factor
    
Make ReadMe with Workflow?
 make_master_file_struct
 preprocess
 process_stats
 - gen_overview_figs
 - gen_unity_plots
 