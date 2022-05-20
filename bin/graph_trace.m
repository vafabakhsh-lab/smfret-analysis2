function graph_trace(time, time_unit, donor, ...
    acceptor, fret, current_trace, file_id)

%{ 

This function ...

    
    Parameters
    ----------
    
    time:
    1-D array of floats corresponding to the time-axis of the trace
    
    time_unit:
    the time resolution of the movie being analyzed and is used for setting
    axis limits
    
    donor:
    N x D array of donor intensity data where N is the number of traces in
    the .traces file and D is the length of the movie

    acceptor:
    N x D array of acceptor intensity data where N is the number of traces
    in the .traces file and D is the length in the movie
    
    fret:
    N x D array of fret values as calculated as acceptor / (donor +
    acceptor)
    
    current_trace:
    intiger index used to iterate over the rows of donor and acceptor for
    calculating fret
    
    file_id:
    integer number corresponding to the film number being analyzed
    
    
    Returns
    -----------
    
    plots of relevant information for the trace
    
    
    Variables
    ----------
    
    pos1/pos2:
    4 element float array dicating the position of the two plots where the
    floats correspond to % of the window being.
     
    
%}

%{ 

    To do:
    
    Try to clean up the plotting section?
    specifically try breaking out the h1 plots into a series of individual
    plots with hold on; and see how it works
    
%}

clf % clear figure to prevent lag caused by plotting

total_intensity = donor + acceptor;

%plot donor, acceptor and fret
pos1 = [0.06, 0.55, 0.68, 0.4];


h1 = subplot('Position', pos1);

title(['File:' file_id '   Molecule:' int2str(current_trace - 1)]);

plot(time,donor(current_trace,:),'g',...
     time,(acceptor(current_trace,:)),'r',...
     time,(total_intensity(current_trace,:) + 800), 'k' );
 
ylim([-25, max(total_intensity(current_trace,:)) + 800]);
xlim([-2, max(time)+100*time_unit]);

grid on;
zoom on;

pos2 = [0.06, 0.1, 0.68, 0.4];
h2 = subplot('Position', pos2);

plot(time, fret(current_trace,:), 'b');

ylim([-0.02, 1.02]);
xlim([-2, max(time)+100*time_unit]);

grid on;
zoom on;
linkaxes([h1,h2], 'x');

end



