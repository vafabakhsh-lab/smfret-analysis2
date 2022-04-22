function [x, donorx, acceptorx] = select_region(trace_length, time_unit)

%{ 

This function ... lets you select a 

    
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
                
disp(newline)
disp('Click twice to select the range.')
[x,~] = ginput(2);

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



