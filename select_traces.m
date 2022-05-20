%{ 

    About

The primary script used for analyzing smFRET data by manually selection
regions of single traces to be saved for further analysis.

    Parameters
    ----------
    
    traces file:
    .traces


    Returns
    ----------
    
    individual files of saved traces based on user input


    Variables
    ----------

    file_directory:
    location of the traces file to be analyzed and is also where the
    individual trace files will be created

    file_id:
    integer number corresponding to the film number being analyzed

    time_unit:
    the time resolution of the movie being analyzed and is used for setting
    axis limits

    particle_positions:

    film_image:
    tiff file that corresponds to the film being analyzed that has been
    previously loaded into matlab
    

    Changelog
    ----------

%}

% Suppress no terminal ; error/warning
%#ok<*NOPTS>

%{ 


Possible things to add:
- Column to keep keep track of the addition or subtraction of background
- gamma correction functionality?
- Add a '(J)ump to particle' option
- switch file format to .csv?

- Rework the structure of the function so that in the loop you 'extract'
the donor, acceptor, and FRET into a current_particle table. This makes it
easier to introduce gamma correction or other functinoality that needs to
work on a specific trace rather than the entire acceptor data

- extract the donor and acceptor (and maybe fret) values for the current
trace into a temporary variable and then many of the downstream functions
wouldn't need to be passed the current_trace variable

- switch from using if-else statements to switch-case statements?


%}


% close any open windows and clear any variables
close all;
fclose('all');
clearvars;


%%%%%%%%%%%%%%%%%%
% Get User Input %
%%%%%%%%%%%%%%%%%%


% input file directory and file name to be analyzed.
file_directory = input('Directory: ', 's');
cd(file_directory);



% get the name of the .traces file to be analyzed
file_id = input('index # of filename [default = 1]  ', 's');
if isempty(file_id)
    file_id = '1';
end
file_name = ['film' file_id '.traces'] 


time_unit = input('Time unit: [default = 0.03 sec] ');
if isempty(time_unit)
    time_unit = 0.03;
end
    
    
leakage = input('Donor leakage correction: [default = 0.085] ');
if isempty(leakage)
  leakage = 0.085;
end
    

bg_subtraction_input = ...
    input('Would you like to subtract the background? y or n [default = n] ', 's');
if bg_subtraction_input == 'y'
    subtract_bg = true;
else
    subtract_bg = false;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load relevant data & files %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load single-molecule data
[donor, acceptor, number_of_traces] = extract_trace_data(file_name);

% correct acceptor for leakage
acceptor = acceptor - (leakage .* donor);

% calculate fret for all traces
%fret = calculate_fret(donor, acceptor, number_of_traces);
fret = acceptor ./ (donor + acceptor);
% correct for rare instance where donor = -acceptor and results in inf
fret(isinf(fret)) = 0;

% get the length of traces
[~, trace_length] = size(donor);

% generate time array
time = (1:trace_length) * time_unit;


% load file image
image_name = ['film' file_id '.tif'] 
film_image = imread(image_name);

% load peak positioins
particle_positions = get_particle_positions(file_id);



%%%%%%%%%%%%%%%%%%
% Begin analysis %
%%%%%%%%%%%%%%%%%%


% iterate over all traces
current_trace = 1;

while current_trace < number_of_traces
    
    % plot trace
    graph_trace(time, time_unit, donor, ...
                acceptor, fret, current_trace, file_id)
            
    hold on;
    show_particle(film_image, particle_positions, current_trace)    
    hold off;
    
    % get user input to advance, go back, or save trace
    option = input('(S)ave, go (B)ack, or skip (ENTER): ', 's');
    
    % validate user input to check that it is empty, s, or b
    while ~isempty(option) && ~contains('sbSB', option)
        disp(['Your input "' option '" is invalid. Please enter a valid option.']);
        option = input([newline ...
            'skip (enter), (s)ave this trace, go (b)ack a trace: '], 's');
    end
     
    
    if isempty(option)
        % go to next trace if no input
        current_trace = current_trace + 1; 
            
    elseif option == 'b' % go back a trace
        % check if at the first trace
        if current_trace - 1 < 1
            disp(newline)
            disp('You are already at the first trace.')
            disp(newline)
        else
            % set index back by 1
            current_trace = current_trace - 1;
        end
    
    
    elseif option == 's'
        % save a trace
        % select part to trace to save
        [x, donorx, acceptorx] = select_region(trace_length, time_unit);
        
        % save the specified region
        save_region(...
            x, donorx, acceptorx, donor, acceptor, ...
            current_trace, file_id, time, subtract_bg ...
            )
        
        % go to next trace
        current_trace = current_trace + 1;
    end
end
