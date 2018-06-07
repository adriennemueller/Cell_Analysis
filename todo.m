%%% FUNCTIONS THAT WANT TO REPLACE process, process2, attIn_att_Out,
%%% aIaOContrasts, WmIn, wmOut etc.

% Will probably have wrapper functions for gen_overviw_fig - to loop
% through all procssd files, and plot dprime unity - with favourite
% selections - or could just put in comment in top of mfile.
function gen_overview_figs( mfs )

end


% Will make and saveout a summary of att-in/att-out, drug-on/drug-off by
% direction for a given file. May include anova sub_fig.
function rslt = gen_overview_fig( file, trial_type )


end




% Will replace original d_prime unitys for vis/combined, etc and plot
% dprimes handed to it.
function rslt = plot_dprime_unity()
end



%%% THINK %%%

% Should combined drug_svm?


% MAKE SVM MORE GENERIC FOR other label types (attin/attout, e.g.)


% Plot SVM performanc as increase number of neurons included

% Plot SVM performance vs scrambled SVM Performance (unity plot? or two
% histograms)