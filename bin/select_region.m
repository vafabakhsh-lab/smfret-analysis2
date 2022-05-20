function [x, donorx, acceptorx] = select_region(trace_length, time_unit)

%{ 

This function facillitates selection a region from a graph that will be
used to save that specific region of a trace as a *.dat file for further 
downstream analysis.

    
    Parameters
    ----------
    
    trace_length:
    the total number of data points in the trace and is used to determine
    if the user clicks outside the bounds of the graph limits

    time_unit:
    time resolution of the film being analyzed and is used to convert the
    x-axis coordinate from time into an integer that can be used as an
    index for saving the selected region
    
    Returns
    -----------
    
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


    Variables
    ----------
    
    temp:
    temporary value for storing one of the elements of x necessary for
    re-ordering the elements of x

    
    
%}
                
disp(newline)
disp('Click twice to select the range.')
[x,~] = ginput(2);

% convert
x(1) = round( x(1) / time_unit );
x(2) = round( x(2) / time_unit );


% Validate user input
if x(1) < 1
    x(1) = 1;
end

if x(2) < 1
    x(2) = 1;
end

if x(1) > trace_length
    x(1) = trace_length;
end

if x(2) > trace_length
    x(2) = trace_length;
end

if x(1) > x(2)
    temp = x(1);
    x(1) = x(2);
    x(2) = temp;
end

disp(['Selected range is frame ', num2str(x(1)), ...
    ' to frame ', num2str(x(2))])
number_of_frames = x(2) - x(1) + 1;


disp('Click twice to select range for calculating donor background.')
[donorx, ~] = ginput(2);
disp('Click twice to select range for calculating acceptor background.')
[acceptorx, ~] = ginput(2);


donorx(1) = round( donorx(1) / time_unit );
donorx(2) = round( donorx(2) / time_unit );

acceptorx(1) = round( acceptorx(1) / time_unit );
acceptorx(2) = round( acceptorx(2) / time_unit);

end



