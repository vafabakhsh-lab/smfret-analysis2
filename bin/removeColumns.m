
clear all;

[file_list, number_of_files] = get_directory();

number_of_columns = input('How many columns are in the source file? ');

destination_folder = input('Where is the data to be saved? ', 's');


for file_index = 1:number_of_files
    
    file_name = file_list(file_index).name;
    
    disp(file_name)
    
    % read data in the trace file
    trace_data = read_trace(file_name, number_of_columns);           
    
    % copy the donor and acceptor data into a new variable
    output = [trace_data.Donor trace_data.Acceptor];
    
    % specify file name the data is to be saved as
    destination_file = strcat(destination_folder,'\', file_list(file_index).name);
    
    % save Donor/Acceptor data
    writematrix(output, destination_file, 'Delimiter', 'comma');     
end

% ========================================================================
% Declare Functions for use in the script
% ========================================================================
function [file_list, number_of_files] = get_directory(~)
    % Get list of *.dat files in a directory
    cd(input('Where are the selected traces? ', 's'));  
    [~, path_name] = uigetfile('*.dat');
    cd(path_name)
    
    file_list = dir('*tr*.dat');
    number_of_files = length(file_list);
end

function [trace_data] = read_trace(file_name, number_of_columns)

    trace_data = readtable(file_name);
    
    % check for the number of columns before naming columns
    
    if number_of_columns == 2
        trace_data.Properties.VariableNames = {'Donor', 'Acceptor'};
    end
    
    if number_of_columns == 4
        trace_data.Properties.VariableNames = {'Time', 'Donor', ...
                                           'Acceptor', 'Gamma'};
    end
    
    if number_of_columns == 6
        trace_data.Properties.VariableNames = {'Time', 'Donor', ...
                                           'Acceptor', 'Gamma', ...
                                           'Donor_BG', 'Acceptor_BG'};
    end
    
end


