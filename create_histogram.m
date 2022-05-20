%{ 

    About

    Parameters & Variables
    ----------------------

    trace_list:
    a list of the trace files in the target directory

    number_of_traces:
    integer number of all files following the *tr*.dat naming format of
    trace files

    bin_size:    
    the size of the bins used when calculating fret histograms
        
    fret_bins:
    1-D array consisting of bin centers
    
    number_bins:
    number of elements in the fret_bins array
    
    max_cc_lag:
    the maximum number of data points away used to calculate the
    correlation between donor and acceptor

    smoothing_window:
    integer number of time points for a moving average used to smooth data

    N_His_f:
    cummulative fret histogram generated from the raw fret trajectory for
    all traces saved

    N_His_fs:
    cummulative fret histogram generated from the smoothed fret trajectory
    for all traces saved

    N_His_F:
    deprecated and no longer used

    N_His_F_gamma:
    ???

    N_His_fret:
    deprecated and no longer used
    
    N_total_XC:
    cummulative cross-correlation data for all traces saved

    N_total_XC_old:
    deprecated and no longer used

    N_total_XC_smooth:
    deprecated and no longer used
    
    skipped_traces:
    cell array to store the name of skipped traces

    skip_count:
    integer count of how many traces were skipped - used for indexing the
    trace name into skipped_traces
    
    particle_count:
    count of how many particles were saved and is used for normalizing the
    data before saving
    
    current_trace:
    integer representing the index of the current trace file being worked
    on and is used to get the specific trace-file name to load

    trace_data:
    multi-dimensional array generated from a single trace file - the number
    of columns can vary but the first three should be time, donor, acceptor
    
    trace_length:
    the number of time points in the current trace_data

    number_of_columns:
    the number of columns in the current trace_data 

    tau:

    
    d_X_a:



    Returns
    -------

    total_fret.csv/.dat:
    csv and dat format of the normalized histogram data for the traces
    saved from a film

        data columns in the file:
            fret: fret bins

            raw: histogram generated from raw fret data

            smoothed: histogram generated from fret data passed through a
            moving average

            filtered-fret: histogram generated from fret data passed
            through a non-linear filter

            filtered-d.a: histogram generated from fret data where the
            donor and acceptor values were passed through a non-linear
            filter
       

    total_xc.csv/.dat:
    csv and dat format of the cross-correlation data for the traces saved
    from a film

        data columns in the file:
            delta-time: the lag time expressed in data-points (aka tau)

            cross-correlation: data from the unfiltered fret normalized by 
            a newer method of (avg_donor + avg_acceptor)



    Changelog
    ----------

%}

% Suppress errors & warnings about:
% no terminal ';' 
% warning about array pre-allocation
% hist is deprecated and use histogram instead

%#ok<*NOPTS>
%#ok<*SAGROW>
%#ok<*HIST>

%{ 
   
Current task: 


Possible things to add or address:
- 'truncation' of the data to be an even multiple of the averaging window
- create chung-kennedy filter
- list of skipped traces for reproducibility
- generate histogram for saved traces and for all traces
- show particle during smoothing?
- the output final should have clear headings
- make sure to retain backwards compatibility with older trace formats
- make the read in better: 
        e.g. trace_data with columns titled time, donor, acceptor


%}


% close any open windows and clear any variables
close all;
fclose('all');
clearvars;


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get list of trace files %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

[FileName,PathName] = uigetfile('*.dat');
cd(PathName);    

trace_list = dir('*tr*.dat')

number_of_traces = length(trace_list)


%%%%%%%%%%%%%%%%%%
% Get User Input %
%%%%%%%%%%%%%%%%%%

timeunit = input('Time unit: [default=0.03 sec] ');
if isempty(timeunit)
    timeunit = 0.03;
end


%%%%%%%%%%%%%%%%%%%%
% define variables %
%%%%%%%%%%%%%%%%%%%%

% histogram bin variables
bin_size = 0.01; 
fret_bins = -0.2:bin_size:1.2;
number_bins = max(size(fret_bins));

max_cc_lag = 60; % max cross-correlation lag (60 datapoints)

smoothing_window = 4; % fret smoothing window

% empty arrays to store cummulutive fret data
N_His_f = zeros(1, number_bins);
N_His_fs = zeros(1, number_bins);
N_His_F = zeros(1, number_bins);
N_His_F_gamma = zeros(1, number_bins);
N_His_fret = zeros(1, number_bins);

% empty arrays to store cummulutive fret data
N_total_XC = zeros(1,max_cc_lag);
N_total_XC_old = zeros(1,max_cc_lag);
N_total_XC_smooth = zeros(1,max_cc_lag);

% misc variables
skipped_traces = {}; % cell array to store what traces were skipped
skip_count = 0; % count of how many traces were skipped
particle_count = 0; % count of how many particles were saved

%%%%%%%%%%%%%%%%%%
% begin analysis %
%%%%%%%%%%%%%%%%%%


for current_trace = 1:number_of_traces

    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % load trace file data %
    %%%%%%%%%%%%%%%%%%%%%%%%
    
    file_name = trace_list(current_trace).name;
    disp(file_name)
    
    import_options = detectImportOptions(file_name);

    trace_data = readtable(file_name, import_options);
    trace_data = table2array(trace_data);
    
    
    [trace_length, number_of_columns] = size(trace_data);
     
    
    % (:,n) = all rows from column n
    time = trace_data(:,1);
    donor = trace_data(:,2);
    acceptor = trace_data(:,3);
    donor_gamma = trace_data(:,4);
    
    % smooth donor and acceptor by averaging over a window
    smooth_donor = smooth(donor, smoothing_window);
    smooth_acceptor = smooth(acceptor, smoothing_window);
    
    % insert non-linear filter function here once created
    % [filtered_donor, filtered_acceptor] = nl_filter(donor, acceptor);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calculate fret and generate histogram for a trace%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  
    fret = acceptor ./ (acceptor + donor);
    smooth_fret = smooth_acceptor ./ (smooth_acceptor + smooth_donor);
    % note: the . before a mathematical operator tells matlab it is an
    % element-by-element calculation and not something else
    


    [n_f, bins] = hist(fret, fret_bins);           
    [n_fs, ~] = hist(smooth_fret, fret_bins); 
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calculate cross-correlation %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    if trace_length > max_cc_lag 
        [d_X_a] = calculate_cross_corr(donor, acceptor);
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plotting & style settings %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    clf; % clear figure
    
    total_intensity = donor + acceptor ;
    smooth_total_intensity = smooth_donor + smooth_acceptor;
    
    % plot donor, acceptor, and total_intensity in window 1
    h1 = subplot(5,2,[1,2,3,4]);
    cla(h1);

    patchline(time,donor, ...
        'edgecolor','g', ...
        'edgealpha',0.4);
    hold on;
    patchline(time, smooth_donor, ...
        'linewidth',1.5, ...
        'edgecolor','g');
    patchline(time, acceptor, ...
        'edgecolor','r', ...
        'edgealpha',0.4);
    patchline(time, smooth_acceptor, ...
        'linewidth',1.5, ...
        'edgecolor','r');
    patchline(time, total_intensity + 150, ...
        'edgecolor','k', ...
        'edgealpha',0.4);
    patchline(time, smooth_total_intensity + 150, ...
        'linewidth',1.5, ...
        'edgecolor','k');
    title(['File: ' file_name '   Molecule:' int2str(current_trace)]);
    ylim([0 20 + max(150 + acceptor + donor)]);
    xlim([min(time) max(time)]);
    grid on; zoom on;

    % plot fret trace in window 2   
    h2 = subplot(5, 2, [5, 6]);
    patchline(time, fret, ...
        'edgecolor', 'c', ...
        'edgealpha', 0.5);
    patchline(time, fret, ...
        'linewidth', 1.5, ...
        'edgecolor', 'b');
    patchline(time, smooth_fret, ...
        'linewidth', 0.4, ...
        'edgecolor', 'b');
    ylim([-0.02 1.02]);
    xlim([min(time) max(time)]);
    grid on; zoom on;
    linkaxes([h1,h2], 'x');

    
    % plot fret histograms in windows 3
    h3 = subplot(5, 2, [7, 9]);
    % get index of most common fret value
    [~, max_fret_index] = max(n_f / sum(n_f)); 
    plot(fret_bins, n_f / sum(n_f), '*'); 
    hold on;
    plot(fret_bins, N_His_f ./ sum(N_His_f), 'r');
    xlim([-0.02 1.02]);
    legend(sprintf('Peak Center = %0.3f', fret_bins(max_fret_index(1))));

    % plot cross-correlation in window 4
    h4 = subplot(5, 2, [8, 10]);
    cla(h4);
    plot(d_X_a_old(1:max_cc_lag), '*b');
    hold on;
    plot(N_total_XC_old / particle_count,'r');
    ylim([-0.6 0.2]);
    legend(sprintf('XC magnitude = %0.3f', d_X_a_old(2)));
      

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % add current trace data to cummulative data or skip %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    option = input('include (enter), skip (s) ', 's');
    if option == 's'
        skip_count = skip_count + 1;
        skipped_traces{skip_count,1} = file_name; 
    else
        particle_count = particle_count + 1;
        
        % normalize the trace histogram by trace length
        n_f  = n_f ./ sum(n_f);
        n_fs = n_fs ./ sum(n_fs);
        
        % add trace-normalized histogram to cummulative histogram
        N_His_f = N_His_f + n_f;
        N_His_fs = N_His_fs + n_fs;
        
        %%%%%%%% cross-correlation  %%%%%%%%%%%%%%%%
        %N_total_XC_old = N_total_XC_old + d_X_a_old(1:max_cc_lag);
        N_total_XC = N_total_XC + d_X_a(1:max_cc_lag);
        %N_total_XC_smooth = N_total_XC_smooth + d_X_a_smooth(1:max_cc_lag);    
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% normalize data by number of particles %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cross-corr normalization
N_total_XC = N_total_XC / particle_count;
N_total_XC_old = N_total_XC_old / particle_count;
N_total_XC_smooth = N_total_XC_smooth / particle_count;

% fret normalization
N_His_f = N_His_f / particle_count;
N_His_fs = N_His_fs / particle_count;
N_His_fret = N_His_fret / particle_count;
N_His_F = N_His_F / particle_count;
N_His_F_gamma = N_His_F_gamma / particle_count;



%%%%%%%%%%%%%%
% save files %
%%%%%%%%%%%%%%

% save cross-corr data
xc_headers = ["delta-time", ...
              "cross-correlation"];

xc_output = [timeunit*(0:max_cc_lag-1)' ...
          N_total_XC'];
      
save('total_XC.dat', 'xc_output', '-ascii') ;

xc_output = array2table(xc_output);
xc_output.Properties.VariableNames = xc_headers;
writetable(xc_output, 'total_XC.csv');


% save fret data
fret_headers = ["fret", ...
                "raw", ...
                "smoothed", ...
                "filtered-fret", ...
                "filtered-d.a", ...
                "gamma"] ;
           
fret_output = [bins' ...
               N_His_f' ...
               N_His_fs' ...
               N_His_fret' ...
               N_His_F' ...
               N_His_F_gamma'];

save('total_FRET.dat', 'fret_output', '-ascii') ;

fret_output = array2table(fret_output);
fret_output.Properties.VariableNames = fret_headers;
writetable(fret_output, 'total_FRET.csv');


% save data on skipped files
skipped_traces = cell2table(skipped_traces);
writetable(skipped_traces, 'skipped-traces.csv')



