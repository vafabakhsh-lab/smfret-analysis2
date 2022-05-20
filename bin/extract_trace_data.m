function [donor, acceptor, number_of_traces] = ...
    extract_trace_data(file_name)

%{ 

This function reads the data from a film.traces file and will return the
donor and acceptor data for all particles in the traces file as well as the
number of particles/traces in the movie.
    
    Parameters
    ----------
    
    file_name:
    file_name of the *.traces file to be analyzed    
    
     
    Returns
    -----------
    
    donor:
    N x D array of donor intensity data where N is the number of traces in
    the .traces file and D is the length of the movie

    acceptor:
    N x D array of acceptor intensity data where N is the number of traces
    in the .traces file and D is the length in the movie
 

    number_of_traces:
    integer number representing the number of trajectories in a .traces
    file 
    
    
    Variables
    ----------
    
    fid:
    file handler matlab uses to read the data
    
    number_of_trajectories:
    combined number of acceptor and donor traces because each trace has two
    trajectories 
    
    number_of_traces:
    number of particles in the *.traces file
    
    raw_data:
    N*2 x D array of data where N is the number of traces, D is the length
    of the movie, the odd rows are donor intensity values, and even rows 
    are acceptor intensity values
    

    Changelog
    -----------
    
    
%}

    
% open file
fid = fopen(file_name,'r');

% get the length of traces
trace_length = fread(fid, 1, 'int32');

% get number of trajectories (both donor and acceptor)
number_of_trajectories = fread(fid, 1, 'int16');

% divide number of trajectories by 2 to obtain the number of
% traces/particles (each trace has 2 trajectories: donor and acceptor)
number_of_traces = number_of_trajectories/2;


% report information to user
disp(newline)
disp(['The number of traces is: ', num2str(number_of_traces)])
disp(['The length of each trace is: ', num2str(trace_length)])
disp(newline)

% read data and close file
raw_data = fread(fid, number_of_trajectories * trace_length, 'int16');
fclose(fid);

% generate empty arrays to store data
donor = zeros(number_of_traces, trace_length);
acceptor = zeros(number_of_traces, trace_length);


% transform raw_data from 1-D array to an N x D array
raw_data = reshape(raw_data, [number_of_trajectories, trace_length]);


for current_trace = 1:(number_of_traces)
    % extract odd rows as donor
    donor(current_trace,:) = raw_data(current_trace*2-1,:);
    % extract even rows as acceptor
    acceptor(current_trace,:) = raw_data(current_trace*2,:); 
end


end

    

% Suppress no terminal ; error/warning
%#ok<*NOPRT>