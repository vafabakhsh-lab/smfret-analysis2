function [fret] = calculate_fret(donor, acceptor, leakage, number_of_traces)

%{ 

This function ...

    
    Parameters
    ----------
    file_name = 
    
    
    Variables
    ----------
    

    Returns
    -----------
    donor
    acceptor
    
    
%}


% correct acceptor  for donor leakage and gammer
corrected_acceptor = gamma.*(acceptor - (leakage * donor));

total_intensity2 = donor + corrected_acceptor;
% make empty array equal to donor and acceptor arrays
fret = zeros(size(donor));

for current_trace = 1:number_of_traces
    
    % calculate total intensity
    total_intensity = donor(current_trace,:) + acceptor(current_trace,:);
    
    % correct acceptor  for donor leakage and gammer
    corrected_acceptor = gamma.*(acceptor(current_trace,:) - ...
                                (leakage * donor(current_trace,:)));
                             
    % calculate fret
    fret(current_trace,:) = corrected_acceptor ./ total_intensity;
    
    % correct any Inf fret values caused by acceptor + donor = 0
    fret(isinf(fret)) = 0;
    
end    
    
    
