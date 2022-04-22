function show_particle(film_image, particle_positions, current_trace)

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

%{ 

The main issue for matlab in reading the .tif saved by smcamera is that
there is no good scaling information available.

To solve this issue, you simply need to tell matlab what the max and min
pixel values will be.

e.g.

imshow(image, [min_intensity max_intensity])

To get those values, you can use:

    [min_intensity, max_intensity] = bounds(image, 'all')

Which will give you the lower and upper bounds for all values in an array
as opposed to each column/row individually.


You can specify specific regions of the image to be displayed by slicing
the image. 

e.g.

imshow(image(y1:y2, x1:x2), [0 max_intensity])

This will show a small section of the image, but it will be small. If you
want the window to be a fixed size then you need to:

% create new figure with specified size 
figure('position', [left, bottom, width, height])  

Or specify image size with:
truesize([100 100])

alternately, can use an 'InitialMagnificiation' parameter!

You can draw a circle using:
drawcircle('Center', [x,y], 'Radius', r);
but it's use is a little hard because it is a movable object and resizable

plot(7,7, 'yo', 'MarkerSize', 50);
%}


% So putting this together you get something like this:
[~, max_intensity] = bounds(film_image, 'all');

% The +1 comes from what I think is the indexing difference in smcamera and
% matlab. where smcamera says the pixel center is 79,259 because it starts
% at 0. So when I write the script to get the co-ordinates from the .pks
% file I need to increment it by one.

donor_x = particle_positions(current_trace,1)+1;
donor_y = particle_positions(current_trace,2)+1;

acceptor_x = particle_positions(current_trace,3)+1;
acceptor_y = particle_positions(current_trace,4)+1;

padding = 6;


% display donor channel of molecule
pos1 = [0.72 0.6 0.3 0.3];
subplot('Position',pos1)
imshow(film_image(donor_y-padding:donor_y+padding, ...
                  donor_x-padding:donor_x+padding...
                  ), ...
       'DisplayRange', [0 max_intensity], ...
       'InitialMagnification', 1500)
title(['Donor (' num2str(donor_x-1) ', ' num2str(donor_y-1) ')']);

% draw bounding circle around particle
viscircles(...
    [7 7], 4,... % set center and radius
    'Color', 'c', ...
    'EnhanceVisibility', false ... % remove white highlight
    );  


% display acceptor channel of molecule
pos2 = [0.72 0.15 0.3 0.3];
subplot('Position', pos2)
imshow(film_image(acceptor_y-padding:acceptor_y+padding, ...
                  acceptor_x-padding:acceptor_x+padding...
                  ), ...
       'DisplayRange', [0 max_intensity], ...
       'InitialMagnification', 1500)
   
% draw bounding circle around particle
viscircles(...
    [7 7], 4,... % set center and radius
    'Color', 'c', ...
    'EnhanceVisibility', false ... % remove white highlight
    );   
   
title(['Acceptor (' num2str(acceptor_x-1) ', ' num2str(acceptor_y-1) ')']);

colormap hot

end
