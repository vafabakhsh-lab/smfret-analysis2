function [] = save_region(...
    x, donorx, acceptorx, donor, acceptor, ...
    current_trace, file_id, time, subtract_bg ...
    )

%{ 

This function saves a trace file based on the a specified set of indices
from the select-region function.

    
    Parameters
    ----------
    
    x:
    two element array of integers that represent x-axis co-ordinates (time
    points) that represent the region of the trace to be saved
    
    donorx:
    two element array of integers that represent x-axis co-ordinates (time
    points) that represent the region of the trace to be used for
    calculating donor background
    
    acceptorx:
    two element array of integers that represent x-axis co-ordinates (time
    points) that represent the region of the trace to be used for
    calculating acceptor background
    
    donor:
    N x D array of donor intensity data where N is the number of traces in
    the .traces file and D is the length of the movie

    acceptor:
    N x D array of acceptor intensity data where N is the number of traces
    in the .traces file and D is the length in the movie
    
    current_trace:
    intiger index used to iterate over the rows of donor and acceptor for
    calculating fret
     
    file_id:
    integer number corresponding to the film number being analyzed
    
    time:
    1-D array of floats corresponding to the time-axis of the trace
    
    subtract_bg:
    logical true or false statement to determine if background is to be
    subtracted from donor and acceptor values

  
    
    Returns
    -----------    
    
    output:
    *.dat file that contains time, donor, acceptor, donor background,
    acceptor background, and background subtraction status for the selected
    region and this file is used to for generating histograms later
    
    
    Variables
    ----------
    
    region_length:
    number of data points in the region
    
    file_name:
    is the name of the file that will be created
    
    mean_donor_background:
    mean value calculated from donorx or entered manually
    
    mean_acceptor_background:
    mean value calculated from acceptorx or entered manually
    
    bg_donor:
    1-D array that containts mean donor background repeated region_length
    times
    
    bg_acceptor:
    1-D array that containts mean acceptor background repeated 
    region_length times
    
    selected_donor:
    donor intensity values from the region to be saved
    
    selected_acceptor:
    donor intensity values from the region to be saved
    
    bg_status:
    numerical true/false logical value that keeps track of whether
    background subtraction has occured
    
    
    
%}


% calculate the length of the region to be saved    
region_length = x(2) - x(1) + 1;
    
% generate file name
file_name = [file_id ' tr' num2str(current_trace - 1) '.dat'];
% current_trace - 1 is to correct for the indexing difference between
% matlab (starts at 1) and smCamera (starts at 0)


% calculate mean donor and acceptor background    

if donorx(1) > donorx(2)
    % get manual input
    mean_donor_background = input('Please enter donor background: ');
else
    mean_donor_background = mean( ...
        donor(current_trace, donorx(1):donorx(2)) );
end

if donorx(1) > donorx(2)
    % get manual input
    mean_acceptor_background = input('Please enter acceptor background: ');
else
    mean_acceptor_background = mean( ...
        acceptor(current_trace, acceptorx(1):acceptorx(2)) );
end


% generate array of background values for saving in the trace file
bg_donor = repmat(mean_donor_background, region_length, 1);
bg_acceptor = repmat(mean_acceptor_background, region_length, 1);


% get donor
selected_donor = donor(current_trace,x(1):x(2));
% get acceptor
selected_acceptor = acceptor(current_trace,x(1):x(2));


if subtract_bg == true
    % get donor
    selected_donor = selected_donor - bg_donor';
    % get acceptor
    selected_acceptor = selected_acceptor - bg_acceptor';
end
    
% generate
bg_status = repmat(subtract_bg, region_length, 1);


output = [...
    time(x(1):x(2))' ...
    selected_donor' ...
    selected_acceptor' ...
    bg_donor ...
    bg_acceptor ...
    bg_status ...
    ];
    
    
save(file_name, 'output', '-ascii');

end



