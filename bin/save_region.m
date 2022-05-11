function [] = save_region(...
    x, donorx, acceptorx, donor, acceptor, ...
    current_trace, file_id, time, subtract_bg ...
    )

%{ 

This function saves a trace file based on the a specified set of indices
from the select-region function.

    
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
%zeros(region_length, 1);



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



