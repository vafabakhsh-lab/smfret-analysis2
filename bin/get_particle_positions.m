function particle_positions = get_particle_positions(file_id)

%{ 

This function will read the *.pks file that has the pixel co-ordinates of
all particles in the smFRET movie and output the [x,y] co-ordinates for
donor and acceptor channels.

    
    Parameters
    ----------
    
    file_id:
    integer number corresponding to the film number being analyzed
    

    Returns
    -----------
    
    particle_positions:
    N x 4 array consisting of the x & y pixel positions of corresponding
    particles
    

    Variables
    ----------

    file_name:
    string value of the .pks file to be loaded

    peak_file:
    string array where each line in the file is a row in the array

    temp_data:
    the .pks file data stripped of the first 6 lines that hold 


    Changelog
    ----------
        
    
%}

% generate file_name based on the file_id being worked on
file_name = ['film' file_id '.pks'];

% read the *.pks file in as a string array
% each row is a line of the file
peak_file = readlines(file_name);

% remove the header lines of the file as well as the empty new line at the
% end of the file
temp_data = peak_file(7:end-1);

% split the rows into their own elements
temp_data = split(temp_data(1:end));

% extract the x,y coordinates for: 
% donor (rows 1 & 3) and acceptor (rows 6 & 8)
particle_positions = [...
                      temp_data(:,1) ...
                      temp_data(:,3) ...
                      temp_data(:,6) ...
                      temp_data(:,8) ...
                     ];

% convert from string to a number (double)                 
particle_positions = str2double(particle_positions);

% round the positions to the nearest pixel
particle_positions = round(particle_positions);

end
             
