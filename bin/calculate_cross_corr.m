function [d_X_a_old, d_X_a, d_X_a_smooth] = ...
    calculate_cross_corr(donor, acceptor, smooth_donor, smooth_acceptor)

%#ok<*AGROW>
%{ 
    About

    Parameters
    ----------

    Returns
    ----------


    Changelog
    ----------

%}

% compute values needed to calculate cross-correlation
avg_donor = mean(donor);
avg_acceptor = mean(acceptor);
avg_smooth_d = mean(smooth_donor); 
avg_smooth_a = mean(smooth_acceptor);

std_donor = std(donor);
std_acceptor = std(acceptor);
std_smooth_d = std(smooth_donor); 
std_smooth_a = std(smooth_acceptor);

delta_donor = donor - avg_donor;
delta_acceptor = acceptor - avg_acceptor;
delta_smooth_d = smooth_donor - avg_smooth_d; 
delta_smooth_a = smooth_acceptor - avg_smooth_a;


trace_length = max(size(donor));

for tau = 0:trace_length 
    sum_tau = 0;
    sum_tau_smooth = 0;
    for t = 1:trace_length-tau
        normalize_number = trace_length-tau;
        sum_tau = ...
            sum_tau + (delta_donor(t) * delta_acceptor(t+tau));
        sum_tau_smooth = ...
            sum_tau_smooth + (delta_smooth_d(t) * delta_smooth_a(t+tau));
    end
    if normalize_number == 0 % why is this here?
        d_X_a(tau+1) = NaN;
    else            
        d_X_a_old(tau+1) = ... 
            (sum_tau / normalize_number) / (std_donor * std_acceptor); 
        % previously it was std / std - check to see if this is correct or
        % not from reza (I think * is the correct way)
        d_X_a(tau+1) = ...
            (sum_tau / normalize_number) / (avg_donor + avg_acceptor);
        d_X_a_smooth(tau+1) = ...
            (sum_tau_smooth / normalize_number)/(avg_smooth_a+avg_smooth_a);
    end
end

% CC equation
% https://wikimedia.org/api/rest_v1/media/math/render/svg/62ffe542d8c069c97cbf975511108207b72c8cd4