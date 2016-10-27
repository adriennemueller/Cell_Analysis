% segregate_stat_struct takes a stat_struct and generates a new stat_struct
% with a subset of the data, depending on what value it is told to filter
% on. eg a value of 'SCH23390' would output a stat_struct with only
% SCH23390 sessions. 

% Valid Drug Values for Filtering
% SKF23390
% SKF81297

% Valid Current Values for Filtering
% 20
% 50 % Make this work for 40 too!
% 100
function rslt = segregate_stat_struct( stat_struct, drug, current )

    rslt = [];

    % Drug Filtering
    if strcmp(drug, 'SCH23390')
        indexes = strcmp([stat_struct.drug], 'SCH23390');
        tmp_rslt = stat_struct(indexes);
    elseif strcmp(drug,'SKF81297')
        indexes = strcmp([stat_struct.drug], 'SKF81297');
        tmp_rslt = stat_struct(indexes);
    else
        disp('Inappropriate value for filtering. Returning original stat_struct.')
        rslt = stat_struct;
    end
    
    % Current Filtering
    for i = 1:length(tmp_rslt)
        if isempty(tmp_rslt(i).attend), continue, end
        currents = [tmp_rslt(i).attend.current];
       
        for j = 1:length(currents)
            if currents(j) == current
                curr_rslt = tmp_rslt(i);
                curr_rslt.attend = tmp_rslt(i).attend(j);
                rslt = [rslt curr_rslt];
            end
        end
    end


end