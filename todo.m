%%% FUNCTIONS THAT WANT TO REPLACE process, process2, attIn_att_Out,
%%% aIaOContrasts, WmIn, wmOut etc.


% Will make and saveout a summary of att-in/att-out, drug-on/drug-off by
% direction for a given file. May include anova sub_fig.
function rslt = gen_overview_fig( file, trial_type )


end


% Will run through all processed files and calculate d's, anovas, etc for
% different epochs of the trial and save them to master_file_struct

%%% SHOULD MAKE PREPROCSS IDENTIFY WHETHER ATTEND OR WM TRIAL, Attend
%%% Contrast, Probe trials
function rslt = process_stats( )

end


% Will replace original d_prime unitys for vis/combined, etc and plot
% dprimes handed to it.
function rslt = plot_dprime_unity()
end



%%% THINK %%%
% Will probably have wrapper functions for gen_overviw_fig - to loop
% through all procssd files, and plot dprime unity - with favourite
% selections - or could just put in comment in top of mfile.

% Should combined drug_svm?


% MAKE SVM MORE GENERIC FOR other label types (attin/attout, e.g.)


% Plot SVM performanc as increase number of neurons included

% Plot SVM performance vs scrambled SVM Performance (unity plot? or two
% histograms)