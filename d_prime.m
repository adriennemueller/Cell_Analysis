% d_prime calculates the d' value for two distributions of data. It takes
% two vectors and calculates the difference in the centers of the
% distributions using the d' metric: (mean1 - mean2) / pooled SD. This
% functions allows for three different ways of calculating the pooled
% standard deviation. All should be valid. It will use the first method by
% default if another method is not specified.
function rslt = d_prime( test_vec, ctrl_vec, sd_type)

    if nargin < 3
        sd_type = 1;
    end

    mu_t = mean(test_vec);
    mu_c = mean(ctrl_vec);
    
    n_t = length(test_vec);
    n_c = length(ctrl_vec);
    
    sd_t = std(test_vec);
    sd_c = std(ctrl_vec);
    
    % All three of these are from the internet and I have no idea which one
    % is best. The first one is from work-learning.com.
    if sd_type == 1
        sd = sqrt( ( (n_t-1)*sd_t^2 + (n_c -1)*sd_c^2) /...
                    (n_t + n_c) );
    elseif sd_type == 2
        sd = sqrt( (sd_t^2 + sd_c^2) / ...
                    2);
    elseif sd_type == 3
        tmp1 = sum(test_vec - mu_t);
        tmp2 = sum(ctrl_vec - mu_c);
        sd = sqrt( tmp1^2 + tmp2^2  / ...
                    n_t + n_c -2 );
    end

    rslt = (mu_t - mu_c) / sd;

end