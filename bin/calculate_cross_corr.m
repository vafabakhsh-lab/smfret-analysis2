function [d_X_a] = calculate_cross_corr(donor, acceptor)

%#ok<*AGROW>
%{ 
    About

    This is an assisting function that calculates the cross-correlation of
    donor and acceptor from a single-molecule.

    Parameters
    ----------
    
    donor: 
    1-d array of N intensity readings, where N is the number of time
    points in the trjacetory

    donor:
    1-d array of N intensity readings, where N is the number of time
    points in the trjacetory

     t:
    current time point (or index) of the 

    tau:
    the lag window for cross correlation calculation between donor(t) and
    acceptor (t+tau)

    sum_tau:
    cummulative value of cross-correlation between donor and acceptor for
    a given tau

    normalize_number:
    integer value that corresponds to the number of ...
    
    trace_length:
    the length of the donor and acceptor arrays


    Returns
    ----------
    
    d_X_a:
    1x60 array that consisting of cross-correlation values for each lag
    window of 1 to 60 time points


    Changelog
    ----------
    

%}

% compute values needed to calculate cross-correlation
avg_donor = mean(donor);
avg_acceptor = mean(acceptor);


std_donor = std(donor);
std_acceptor = std(acceptor);


delta_donor = donor - avg_donor;
delta_acceptor = acceptor - avg_acceptor;


trace_length = max(size(donor));

for tau = 0:trace_length 
    sum_tau = 0;
    for t = 1:trace_length-tau
        normalize_number = trace_length-tau;
        sum_tau = ...
            sum_tau + (delta_donor(t) * delta_acceptor(t+tau));

    end
    if normalize_number == 0 % why is this here?
        d_X_a(tau+1) = NaN;
    else            
        d_X_a(tau+1) = ... 
            (sum_tau / normalize_number) / (std_donor * std_acceptor); 

    end
end

% CC equation
% https://wikimedia.org/api/rest_v1/media/math/render/svg/62ffe542d8c069c97cbf975511108207b72c8cd4