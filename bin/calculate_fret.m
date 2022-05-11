function [fret] = calculate_fret(donor, acceptor, number_of_traces)

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



% make empty array equal to donor and acceptor arrays
fret = zeros(size(donor));

for current_trace = 1:number_of_traces
    
    % calculate total intensity
    total_intensity = donor(current_trace,:) + acceptor(current_trace,:);
    
    % calculate fret
    fret(current_trace,:) = acceptor(current_trace,:) ./ total_intensity;
    
    % correct any Inf fret values caused by acceptor + donor = 0
    fret(isinf(fret)) = 0;
    
end

end
    
    
