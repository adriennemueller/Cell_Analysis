% get_window identifies the event codes that delimit the specified window.
% A correct trial from the data_struct is passed in, as well as string
% identifying what window is desired:
% 'attend'      - cue onset to blank onset
% 'attContrast' - cue onset to blank onset
% 'wm'          - target offset to fixation offset % VERIFY
% 'visual'      - targets onset to either cue onset (attend trials) or target offset (wm trials)
% 'fixation'    - fixation onset to targets onset 
% 'fullNoMotor' - fixation onset to blank onset (attend trials) or to fixation offset (wm trials)
% 
function window = get_window( correct_trial, window_string )
    attend_earlysession_flag = find(correct_trial.event_codes == 121);
    attend_latesession_flag  = find(correct_trial.event_codes == 133);
    wm_flag                  = find(correct_trial.event_codes == 153);

    if strcmp( window_string, 'attend' ) || strcmp( window_string, 'attContrast' )
        if attend_earlysession_flag, window = [121 126]; % Need a function for this because I changed the event codes in Jan/Mar 2016 
        else, window = [133 126];
        end
    elseif strcmp( window_string, 'wm' )
        window = [155 161]; 
    elseif strcmp( window_string, 'visual' )
        if attend_earlysession_flag, window = [124 121]; % Attend Trials Early Sessions
        elseif attend_latesession_flag, window = [124 133]; % Attend Trials Late Sessions
        else, window = [153 155]; % WM Trials
        end
    elseif strcmp( window_string, 'fixation' )
        if wm_flag, window = [120 153]; % WM Trials
        else, window = [120 124]; % Attend Trials
        end
    elseif strcmp( window_string, 'fullNoMotor' )
        if wm_flag, window = [120 161];
        else, window = [120 126]; 
        end
    end
end