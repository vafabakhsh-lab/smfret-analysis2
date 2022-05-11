%{ 

    About

    Parameters
    ----------


    Returns
    ----------


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

Pseudocoding

1. clean workspace
2. get user input
3. get list of trace files
4. declare variables
5. for each trace file in folder
    a. calculate FRET
    b. calculate trace histogram
    c. calculate XC
    d. graph data
    e. save or skip? (no going back)
    f. add trace to sum of histograms
    g. increment number of traces added (for normalization)
6. save histogram data


    
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

trace_list=dir('*tr*.dat')

number_of_traces = length(trace_list)


%%%%%%%%%%%%%%%%%%
% Get User Input %
%%%%%%%%%%%%%%%%%%

timeunit = input('Time unit: [default=0.03 sec] ');
if isempty(timeunit)
    timeunit=0.03;
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

    trace_list(current_trace).name;
    
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
    
    % calculate fret
    fret = acceptor ./ (acceptor+donor);
    smooth_fret = smooth_acceptor ./ (smooth_acceptor + smooth_donor);

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % generate FRET histogram %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%

    [n_f, bins] = hist(fret, fret_bins);           
    [n_fs,~] = hist(smooth_fret, fret_bins); 
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calculate cross-correlation %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    if trace_length > max_cc_lag %
        [d_X_a_old, d_X_a, d_X_a_smooth] = ...
            calculate_cross_corr( donor, acceptor, ...
                                  smooth_donor, smooth_acceptor );
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
    patchline(time, ...
        total_intensity + 150, ...
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
    plot(fret_bins,n_f / sum(n_f),'*'); 
    hold on;
    plot(fret_bins, N_His_f ./ sum(N_His_f), 'r');
    xlim([-0.02 1.02]);
    legend(sprintf('Peak Center = %0.3f', fret_bins(max_fret_index(1))));

    % plot cross-correlation in window 5
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
        N_total_XC_old = N_total_XC_old + d_X_a_old(1:max_cc_lag);
        N_total_XC = N_total_XC + d_X_a(1:max_cc_lag);
        N_total_XC_smooth = N_total_XC_smooth + d_X_a_smooth(1:max_cc_lag);    
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
    "unfiltered-data", ...
    "old-normalization", ...
    "filtered-data"];

xc_output = [timeunit*(0:max_cc_lag-1)' ...
          N_total_XC' ...
          N_total_XC_old' ...
          N_total_XC_smooth'];
      
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



