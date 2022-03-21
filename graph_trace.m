function graph_trace(time, time_unit, donor, ...
    acceptor, fret, current_trace, file_id)

%{ 

This function ...

    
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



clf % clear figure to prevent lag caused by plotting

total_intensity = donor + acceptor;

%plot donor, acceptor and fret
pos1 = [0.08 0.55 0.6 0.4];
h1 = subplot('Position',pos1);

plot(time,donor(current_trace,:),'g',...
     time,(acceptor(current_trace,:)),'r',...
     time,(total_intensity(current_trace,:) + 800), 'k' );
 
ylim([-25, max(total_intensity(current_trace,:)) + 400]);

xlim([-2, max(time)+100*time_unit]);

title(['File:' file_id '   Molecule:' int2str(current_trace-1)]);
grid on;
zoom on;


pos2 = [0.08 0.1 0.6 0.4];
h2 = subplot('Position', pos2);
plot(time,fret(current_trace,:),'b');
ylim([-0.02, 1.02]);
xlim([-2, max(time)+100*time_unit]);

grid on;
zoom on;
linkaxes([h1,h2], 'x');

end



