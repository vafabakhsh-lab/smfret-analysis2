function show_particle(film_image, particle_positions, current_trace)

%{ 

This function will display a specific region of an image to show the
current particle that is currently being worked on.

    
    Parameters
    ----------
    
    film_image:
    tiff file that corresponds to the film being analyzed that has been
    previously loaded into matlab

    particle_positions:
    4-element array containing the x & y coordinates for both donor and
    acceptor channel of a particle

    current_trace:
    intiger index used to iterate over the rows of donor and acceptor for
    calculating fret
        

    Returns
    -----------
    two plots showing the particle currently being analyzed in the donor
    and acceptor channel
    

    Variables
    ----------
    
    max_intensity:
    maximum intensity in the tiff that

    particle_size:
    the radius of the particles - typically what is used in smCamera during
    analysis and is ~4 pixels

    circle_thickness:
    determines the thickness of the circle surrounding the particle    

    padding:
    integer value used to create the region around the particle to be
    displayed

    pos1/pos2:
    4 element float array dicating the position of the two plots where the
    floats correspond to % of the window being. It follows the format of 
    [left bottom width height] where the first two values determine the
    bottom left corner of the figure.

    

    
    
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
[min_intensity, max_intensity] = bounds(film_image, 'all');

% The +1 comes from what I think is the indexing difference in smcamera and
% matlab. where smcamera says the pixel center is 79,259 because it starts
% at 0. So when I write the script to get the co-ordinates from the .pks
% file I need to increment it by one.

donor_x = particle_positions(current_trace, 1) + 1;
donor_y = particle_positions(current_trace, 2) + 1;

acceptor_x = particle_positions(current_trace, 3) + 1;
acceptor_y = particle_positions(current_trace, 4) + 1;

particle_size = 4;
circle_thickness = 2.5;
padding = 6;


% display donor channel of molecule
pos1 = [0.74 0.6 0.3 0.3];
f1 = subplot('Position', pos1);
imshow(film_image(donor_y - padding:donor_y + padding, ...
                  donor_x - padding:donor_x + padding...
                  ), ...
       'DisplayRange', [min_intensity, max_intensity], ...
       'InitialMagnification', 1500, ...
       'Border', 'tight');
title(['Donor (' num2str(donor_x - 1) ', ' num2str(donor_y - 1) ')']);

% draw bounding circle around particle
viscircles(...
    [padding + 1, padding + 1], particle_size,... % set center and radius
    'Color', 'y', ...
    'LineWidth', circle_thickness, ...
    'EnhanceVisibility', false ... % remove white highlight
    );  


% display acceptor channel of molecule
pos2 = [0.74 0.15 0.3 0.3];
f2 = subplot('Position', pos2);
imshow(film_image(acceptor_y - padding:acceptor_y + padding, ...
                  acceptor_x - padding:acceptor_x + padding...
                  ), ...
       'DisplayRange', [min_intensity, max_intensity], ...
       'InitialMagnification', 1500, ...
       'Border', 'tight');
   
% draw bounding circle around particle
viscircles(...
    [padding + 1, padding + 1], particle_size,... 
    'Color', 'y', ...
    'LineWidth', circle_thickness, ...
    'EnhanceVisibility', false ... 
    );   
   
title(['Acceptor (' num2str(acceptor_x - 1) ', ' num2str(acceptor_y - 1) ')']);

linkaxes([f1, f2], 'xy');


% specify colormap based on smCamera software colormap

map =   [0         0         0
         0         0    0.0196
         0         0    0.0392
         0         0    0.0588
         0         0    0.0784
         0         0    0.1020
         0         0    0.1216
         0         0    0.1412
         0         0    0.1608
         0         0    0.1804
         0         0    0.2039
         0         0    0.2235
         0         0    0.2627
         0         0    0.2824
         0         0    0.3059
         0         0    0.3255
         0         0    0.3451
         0         0    0.3647
         0         0    0.3843
         0         0    0.4078
         0         0    0.4275
         0         0    0.4471
         0         0    0.4667
         0         0    0.4863
         0         0    0.5294
         0         0    0.5490
         0         0    0.5686
         0         0    0.5882
         0         0    0.6118
         0         0    0.6314
         0         0    0.6510
         0         0    0.6706
         0         0    0.6902
         0         0    0.7137
         0         0    0.7333
         0         0    0.7529
         0         0    0.7922
         0         0    0.8157
         0         0    0.8353
         0         0    0.8549
         0         0    0.8745
         0         0    0.8941
         0         0    0.9176
         0         0    0.9373
         0         0    0.9569
    0.0157         0    0.9765
    0.0353         0    1.0000
    0.0745         0    0.9608
    0.0902         0    0.9373
    0.1098         0    0.9176
    0.1294         0    0.8941
    0.1490         0    0.8745
    0.1647         0    0.8549
    0.1843         0    0.8314
    0.2039         0    0.8118
    0.2235         0    0.7882
    0.2392         0    0.7686
    0.2588         0    0.7451
    0.2784         0    0.7255
    0.3176         0    0.6824
    0.3176         0    0.6627
    0.3176         0    0.6392
    0.3176         0    0.6196
    0.3176         0    0.5961
    0.3176         0    0.5765
    0.3176         0    0.5569
    0.3176         0    0.5333
    0.3137         0    0.5137
    0.3137         0    0.4902
    0.3137         0    0.4706
    0.3137         0    0.4471
    0.3137         0    0.4078
    0.3137         0    0.3843
    0.3098         0    0.3647
    0.3294         0    0.3412
    0.3490         0    0.3216
    0.3686         0    0.2980
    0.3882         0    0.2784
    0.4078         0    0.2588
    0.4275         0    0.2353
    0.4471         0    0.2157
    0.4667         0    0.1922
    0.4863         0    0.1725
    0.5255         0    0.1294
    0.5451         0    0.1098
    0.5647         0    0.0863
    0.5843         0    0.0667
    0.6039         0    0.0431
    0.6235         0    0.0235
    0.6431         0         0
    0.6627         0         0
    0.6824         0         0
    0.7059         0         0
    0.7255         0         0
    0.7686         0         0
    0.7882         0         0
    0.8078         0         0
    0.8314         0         0
    0.8510         0         0
    0.8706         0         0
    0.8941         0         0
    0.9137         0         0
    1.0000         0         0
    1.0000         0         0
    1.0000         0         0
    1.0000         0         0
    1.0000    0.0392         0
    1.0000    0.0627         0
    1.0000    0.0824         0
    1.0000    0.1059         0
    1.0000    0.1255         0
    1.0000    0.1451         0
    1.0000    0.1686         0
    1.0000    0.1882         0
    1.0000    0.2118         0
    1.0000    0.2314         0
    1.0000    0.2510         0
    1.0000    0.2745         0
    1.0000    0.3176         0
    1.0000    0.3333    0.0157
    1.0000    0.3529    0.0353
    1.0000    0.3725    0.0549
    1.0000    0.3922    0.0745
    1.0000    0.4118    0.0941
    1.0000    0.4275    0.1098
    1.0000    0.4471    0.1294
    1.0000    0.4667    0.1490
    1.0000    0.4863    0.1686
    1.0000    0.5059    0.1882
    1.0000    0.5255    0.2078
    1.0000    0.5608    0.2431
    1.0000    0.5804    0.2627
    1.0000    0.6000    0.2824
    1.0000    0.6196    0.3020
    1.0000    0.6392    0.3216
    1.0000    0.6392    0.3020
    1.0000    0.6392    0.2784
    1.0000    0.6392    0.2549
    1.0000    0.6392    0.2314
    1.0000    0.6392    0.2078
    1.0000    0.6392    0.1843
    1.0000    0.6392    0.1412
    1.0000    0.6392    0.1176
    1.0000    0.6392    0.0941
    1.0000    0.6392    0.0706
    1.0000    0.6392    0.0471
    1.0000    0.6392    0.0235
    1.0000    0.6392         0
    1.0000    0.6392         0
    1.0000    0.6392         0
    1.0000    0.6392         0
    0.9725    0.6392         0
    0.9412    0.6392         0
    0.8824    0.6392         0
    0.8510    0.6392         0
    0.8196    0.6392         0
    0.7922    0.6392         0
    0.7608    0.6392         0
    0.7294    0.6392         0
    0.7020    0.6392         0
    0.6706    0.6392         0
    0.6392    0.6392         0
    0.6588    0.6392         0
    0.6784    0.6392         0
    0.6980    0.6627    0.0118
    0.7373    0.7098    0.0353
    0.7569    0.7333    0.0471
    0.7765    0.7569    0.0627
    0.7961    0.7804    0.0745
    0.8196    0.8039    0.0863
    0.8392    0.8314    0.0980
    0.8588    0.8549    0.1137
    0.8784    0.8784    0.1255
    0.8980    0.9020    0.1373
    0.9176    0.9255    0.1490
    0.9373    0.9490    0.1608
    0.9569    0.9725    0.1765
    1.0000    1.0000    0.2000
    1.0000    1.0000    0.2118
    1.0000    1.0000    0.2275
    1.0000    1.0000    0.2392
    1.0000    1.0000    0.2510
    1.0000    1.0000    0.2627
    1.0000    1.0000    0.2784
    1.0000    1.0000    0.2902
    1.0000    1.0000    0.3020
    1.0000    1.0000    0.3137
    1.0000    1.0000    0.3255
    1.0000    1.0000    0.3529
    1.0000    1.0000    0.3647
    1.0000    1.0000    0.3765
    1.0000    1.0000    0.3922
    1.0000    1.0000    0.4039
    1.0000    1.0000    0.4157
    1.0000    1.0000    0.4275
    1.0000    1.0000    0.4392
    1.0000    1.0000    0.4549
    1.0000    1.0000    0.4667
    1.0000    1.0000    0.4784
    1.0000    1.0000    0.4902
    1.0000    1.0000    0.5176
    1.0000    1.0000    0.5294
    1.0000    1.0000    0.5412
    1.0000    1.0000    0.5569
    1.0000    1.0000    0.5686
    1.0000    1.0000    0.5804
    1.0000    1.0000    0.5922
    1.0000    1.0000    0.6039
    1.0000    1.0000    0.6196
    1.0000    1.0000    0.6314
    1.0000    1.0000    0.6431
    1.0000    1.0000    0.6549
    1.0000    1.0000    0.6824
    1.0000    1.0000    0.6941
    1.0000    1.0000    0.7059
    1.0000    1.0000    0.7176
    1.0000    1.0000    0.7333
    1.0000    1.0000    0.7451
    1.0000    1.0000    0.7569
    1.0000    1.0000    0.7686
    1.0000    1.0000    0.7843
    1.0000    1.0000    0.7961
    1.0000    1.0000    0.8078
    1.0000    1.0000    0.8196
    1.0000    1.0000    0.8471
    1.0000    1.0000    0.8588
    1.0000    1.0000    0.8706
    1.0000    1.0000    0.8824
    1.0000    1.0000    0.8980
    1.0000    1.0000    0.9098
    1.0000    1.0000    0.9216
    1.0000    1.0000    0.9333
    1.0000    1.0000    0.9490
    1.0000    1.0000    0.9608
    1.0000    1.0000    0.9725
    1.0000    1.0000    1.0000];

% apply colormap
colormap(map)


end
