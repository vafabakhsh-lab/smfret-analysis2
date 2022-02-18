%{ 

    About

    Parameters
    ----------

    Returns
    ----------


    Changelog
    ----------

%}
% Suppress no terminal ; error/warning
%#ok<*NOPTS>

%{ 

    Pseudocoding
    
    1. clean workspace
    2. declare variables
    3. get data directory
    4. get data parameters
    5. read data
    6. for each trace in data 
        a. calculate FRET
       > b. graph trace
        c. ask to save or skip
    7. save trace
        a. select part of trace to save
        b. define background
        c. write data to file
    

%}


% close any open windows and clear any variables
close all;
fclose('all');
clearvars;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get User Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% input file directory and file name to be analyzed.
pth=input('Directory: ', 's');
cd(pth);

% Specify how many points to average.
sm=3;
% Specify the Gamma factor
gamma=1.0;

gap=0:0.01:1;
sg=size(gap);sg=sg(2);
N_His=zeros(1,sg);
N_His_raw=zeros(1,sg);

% get the traces file name to be analyzed
file_id = input('index # of filename [default=1]  ', 's');
if isempty(file_id)
    file_id='1';
end
file_name = ['film' file_id '.traces'] 


time_unit=input('Time unit: [default=0.06 sec] ');
if isempty(time_unit)
    time_unit=0.06;
end
    
    
leakage=input('Donor leakage correction: [default=0.09] ');
if isempty(leakage)
  leakage=0.09;
end
    
    
bg_definition_input=input('Manually define background? y or n [default=y] ');
if isempty(bg_definition_input)
    bg_definition_input='y';
end
    
bg_definition=1;

if bg_definition_input=='n'
    bg_definition_input=0;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[donor, acceptor, number_of_traces] = extract_trace_data(file_name);

[~, trace_length] = size(donor)

time = (1:trace_length) * time_unit;



fret = calculate_fret(donor, acceptor, leakage, number_of_traces)


while current_trace < number_of_traces
    
    current_trace = current_trace + 1;
    
    graph_trace(time, time_unit, donor, acceptor, fret, current_trace)
    
    % advance, go back, or save trace
    option = input('skip (enter), (s)ave this trace, go (b)ack a trace: ', 's');
    % validate user input
    while ~isempty(option) && ~contains('sb', option)
        disp(['Your input "' option '" is invalid. Please enter a valid option.']);
        option = input([newline ...
            'skip (enter), (s)ave this trace, go (b)ack a trace: '], 's');
    end
    
    if option == 'b'
        % set index back by 2 - to account for the +1 in the loop
        current_trace = current_trace - 2;
    end
    
    if option =='s'
        % select part to trace
        % select_region()
    end
    

end
