function [fret] = calculate_fret(donor, acceptor, number_of_traces)

%{ 

This function was created to take the bulk donor and acceptor data from all
traces in a .traces file that has been read into matlab and return the fret
values for each particles' trajectory.

    
    Parameters
    ----------

    donor:
    N x D array of donor intensity data where N is the number of traces in
    the .traces file and D is the length of the movie

    acceptor:
    N x D array of acceptor intensity data where N is the number of traces
    in the .traces file and D is the length in the movie
 

    number_of_traces:
    integer number representing the number of trajectories in a .traces
    file 
    

    Returns
    -----------
    
    fret:
    N x D array of fret values as calculated as acceptor / (donor +
    acceptor)


    Variables
    ----------
    
    current_trace:
    intiger index used to iterate over the rows of donor and acceptor for
    calculating fret
    
    total_intensity:
    donor + acceptor intensities used to make fret calculation more clear   
    

    Changelog
    -----------  
  

%}



% make empty array equal to donor and acceptor arrays
fret = zeros(size(donor));

for current_trace = 1:number_of_traces
    
    % calculate total intensity
    total_intensity = donor(current_trace, :) + acceptor(current_trace,:);
    
    % calculate fret
    fret(current_trace, :) = acceptor(current_trace, :) ./ total_intensity;
    
    % correct any Inf fret values caused by acceptor + donor = 0
    fret(isinf(fret)) = 0;
    
end

end
    
    
