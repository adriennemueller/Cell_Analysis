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
% Will return a pair of numbers:
% end_code of the window
% win_length - how many samples should be counted back from that window

%%% WM CODES:
% 9 / 9 / 9 Trial Start
% 120 Fix Start
% 153 Cue On
% 155 Cue Off / Delay Start
% 157 Fix Off
% 161 Acquire Target/Reward
% 18 / 18 / 18 Trial End
function win_info = get_window( e_codes, e_times, window_string )
    attend_earlysession_flag = find(e_codes == 121);
    attend_latesession_flag  = find(e_codes == 133);
    wm_flag                  = find(e_codes == 153);

    if ismember( window_string, {'attend', 'attContrast'} )
        if attend_earlysession_flag, trial_window = [121 126]; % Need a function for this because I changed the event codes in Jan/Mar 2016 
        else, trial_window = [133 126];
        end
    elseif strcmp( window_string, 'blank' )
        trial_window = [126 128];
    elseif strcmp( window_string, 'post_blank' )
        trial_window = [128 128];    
    elseif strcmp( window_string, 'wm_delay' )
        trial_window = [155 157];
    elseif strcmp (window_string, 'wm_response' )% These are variable durations.
         trial_window = [157 161];
    elseif ismember( window_string, {'attend_visual', 'wm_visual', 'visual'} )
        if attend_earlysession_flag, trial_window = [124 121]; % Attend Trials Early Sessions
        elseif attend_latesession_flag, trial_window = [124 133]; % Attend Trials Late Sessions
        else, trial_window = [153 155]; % WM Trials
        end
    elseif ismember( window_string, {'attend_fixation', 'wm_fixation', 'fixation'} )
        if wm_flag, trial_window = [120 153]; % WM Trials
        else, trial_window = [120 124]; % Attend Trials
        end
    elseif ismember( window_string, {'attend_fullNoMotor', 'wm_fullNoMotor','fullNoMotor'} ) % Up to and including Blank
        if wm_flag, trial_window = [120 161];
        else, trial_window = [120 128]; %[120 126]; 
        end
    elseif ismember( window_string, {'attend_reward', 'wm_reward','reward'} )
        if wm_flag, trial_window = [157 161];% WM Trials
        else, trial_window = [132 132]; % No actual end for this;
        end
    end
    
    end_code = trial_window(2);
    win_length = get_win_length( e_codes, e_times, trial_window, window_string );
    
    win_info = [end_code win_length];
end


function win_length = get_win_length( e_codes, e_times, trial_window, window_string )

    %attend_earlysession_flag = find(correct_trial.event_codes == 121);
    %attend_latesession_flag  = find(correct_trial.event_codes == 133);
    wm_flag                  = find(e_codes == 153);

    fix_win_length = 300;  % Always take the last 300ms of the fixation window
    wm_resp_win_length  = 300;  % Always take the last 500ms of the delay window in the WM task
    reward_length  = -300; % Take the 300ms going FORWARD from the end of the trial
    post_blank_length = -300;
    
    % If fixation window 
    if ismember(window_string, {'fixation', 'attend_fixation', 'wm_fixation'} )
        win_length = fix_win_length; 
        
    % Useful because these ranges are variable duration.    
    elseif strcmp( window_string, 'wm_response' ) 
        win_length = wm_resp_win_length;
        
    elseif strcmp( window_string, 'reward' )
        win_length = reward_length;

    elseif strcmp( window_string, 'post_blank' )
        win_length = post_blank_length;

    % fullNoMotor window
    elseif strcmp( window_string, 'fullNoMotor' )
        if wm_flag, fix_end_code = 153;
        else, fix_end_code = 124; 
            end_time = e_times( e_codes == trial_window(2) );
            fix_end_time = e_times( e_codes == fix_end_code );
            win_length = (end_time - fix_end_time) + fix_win_length;
        end
        
    % If visual or attend or working memory windows - consistent times
    else
        beg_time = e_times( e_codes == trial_window(1) );
        end_time = e_times( e_codes == trial_window(2) );
        win_length = end_time - beg_time;
    end
    
end