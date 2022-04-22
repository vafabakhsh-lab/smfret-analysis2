
clear;
close all;

path=input('where are the selected traces?  ', 's');
cd(path);
[FileName,PathName] = uigetfile('*.dat');
cd(PathName);


      
A=dir('*tr*.dat');
numberfiles=length(A);

max_length = 0;

for file_index = 1:numberfiles
    
    A(file_index).name;
    fname=A(file_index).name(1:end-4);
    fid=fopen(A(file_index).name,'r');
    TraceFile = fscanf(fid,'%g %g %g %g %g %g',[6 inf]);
    fclose(fid);
    trace_size = size(TraceFile);
    trace_length = trace_size(2);
    
    if trace_length > max_length
        
        max_length = trace_length;
        
    end 
end

fret_traces = NaN(max_length, numberfiles);
fret_and_time = NaN(max_length + 1, numberfiles);

for file_index = 1:numberfiles
   
    A(file_index).name;
    fname=A(file_index).name(1:end-4);
    fid=fopen(A(file_index).name,'r');
    TraceFile = fscanf(fid,'%g %g %g %g %g %g',[6 inf]);
    fclose(fid);
    trace_size = size(TraceFile);
    trace_length = trace_size(2);

    donor=+TraceFile(2,:)';       
    acceptor=+TraceFile(3,:)';
    fret=acceptor./(donor+acceptor); 

    fret_traces(1:trace_length, file_index) = ...
        fret;

    
    fret_and_time(1:trace_length + 1, file_index) = ...
        [trace_length; fret];

end



%mymap = load('D:\smFret\programs2\SchamberPrograms\FretHeatmapColor.mat')

fret_traces = fret_traces'; % Transposes RawFret
fret_and_time = fret_and_time';

meanfret = mean(fret_traces, 2, 'omitnan'); % calculates mean of RawFret by row (2)
fret_b = [fret_traces meanfret]; % appends MeanFret to the end of RawFret
[~,idx1] = sort(fret_b(:,end)); % sorts the ids of F2 by the last column (MeanFret)
fret_sorted = fret_b(idx1,:); % Rearrange F2 by the sorted ids (idx1)


[~,idx1] = sort(fret_and_time(1:end,1));
fret_time_sorted = fret_and_time(idx1,:);
fret_time_sorted = fret_time_sorted(1:end,2:end);


pos1 = [0.08 0.55 0.9 0.4];
h1 = subplot('Position',pos1);

imagesc(fret_sorted);
set(gca,'xtick',[]);
caxis([-0.1 1]);
colorbar;
title('Traces sorted by mean FRET');

pos2 = [0.08 0.1 0.9 0.4];
h2 = subplot('Position', pos2);

imagesc(fret_time_sorted); 
caxis([-0.1 1]);
colorbar;
title('Traces sorted by length');

linkaxes([h1,h2], 'x');
colormap([0.1 0.1 0.1; turbo(256)])

exportgraphics(gcf, 'trace_heatmap.png', 'Resolution', 1200)