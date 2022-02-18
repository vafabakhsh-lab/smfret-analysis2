function graph_trace(time, time_unit, donor, acceptor, fret, current_trace)

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
h1=subplot(5,1,[1,2,3]);
plot(time,donor(current_trace,:),'g',...
     time,(acceptor(current_trace,:)),'r',...
     time,(total_intensity(current_trace,:) + 800), 'k' );
ylim([-25 max((acceptor(current_trace,:)+donor(current_trace,:)+400))]);
xlim([-20*time_unit max(time)+100*time_unit]);

title(['File:' 1 '   Molecule:' int2str(current_trace-1)]);
grid on;
zoom on;
h2=subplot(2,1,2);
plot(time,fret(current_trace,:),'b');
 ylim([-0.02 1.02]);
 xlim([-20*time_unit max(time)+100*time_unit]);

grid on;
zoom on;
linkaxes([h1,h2], 'x');

