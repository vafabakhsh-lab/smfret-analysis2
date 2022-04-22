function filtered_rawTrace = ckFilter(rawTrace, windows, M, p)

%{ 

This function returns data that has been filtered with a chunk-kennedy
nonlinear filter. Source: https://doi.org/10.1016/0165-0270(91)90118-J

    Parameters
    ----------
    rawTrace: numpy 1-D array
        rawTraceectory
    windows: list of ints
        lengths of windows to calculate forward and reverse predictors
    M: int
      size of window over which the predictors are to be compared
    p: int
       weighting factor
    Returns
    -----------
    filtered_rawTrace: numpy 1-D arraay
        NLF filtered rawTraceectory
%}

time = length(rawTrace);
%
%Set algorithm parameters (discussion on this in Bioinformatics paper - Supplement A):
%
% ^^^ NOTE: User can change following parameter to improve fit ^^^
%
windows = [4 8 16 32]; %Sampling window length
NumberWindows = length(windows); %Number of sampling windows
M = 10; %Weighting parameter length
P = 40; %Parameter that effects sharpness of steps
%
% ——- The Forward-Backward non-linear Algorithm ———-
%
%First - Run each time step and store the average forward and backward
%predictors for each of the sampling windows
%
I_avg_f = zeros(time,NumberWindows);
I_avg_b = zeros(time,NumberWindows);
for g = 1:time
%Run for each sampling window
    for k = 1:NumberWindows
    %Solve for the average forward predictor (w/ some
    %logic to help solve the beginning of the trace
    window = windows(k);
        if g == 1
            I_avg_f(g,k) = rawTrace(1,1);
        elseif g - window - 1 < 0
            I_avg_f(g,k) = sum(rawTrace(1:g-1,1))/g;
        else
            epoint = g - window;
            spoint = g - 1;
            I_avg_f(g,k) = sum(rawTrace(epoint:spoint,1))/window;
        end
%Now do the same for the backward predictor
        if g == time
            I_avg_b(g,k) = rawTrace(g,1);
        elseif g + window > time
            sw = time - g;
            I_avg_b(g,k) = sum(rawTrace(g+1:time,1))/sw;
        else
            epoint = g + window;
            spoint = g + 1;
            I_avg_b(g,k) = sum(rawTrace(spoint:epoint,1))/window;
        end
    end
end
%
%Second - Solve for non-normalized forward and backward weights:
f = zeros(time,NumberWindows);
b = zeros(time,NumberWindows);
for i = 1:time
    for k = 1:NumberWindows
        Mstore_f = zeros(M,1);
        Mstore_b = zeros(M,1);
        for j = 0:M-1
            t_f = i - j;
            t_b = i + j;
            if t_f < 1
                Mstore_f(j+1,1) = (rawTrace(i,1) - I_avg_f(i,k))^2;
            else
                Mstore_f(j+1,1) = (rawTrace(t_f,1) - I_avg_f(t_f,k))^2;
            end
            if t_b > time
                Mstore_b(j+1,1) = (rawTrace(i,1) - I_avg_b(i,k))^2;
            else
                Mstore_b(j+1,1) = (rawTrace(t_b,1) - I_avg_b(t_b,k))^2;
            end
        end
        f(i,k) = sum(Mstore_f)^(-P);
        b(i,k) = sum(Mstore_b)^(-P);
    end
end
% Third - Solve for vector of normalization factors for the weights
C = zeros(time,1);
for i = 1:time
    NumberWindowsstore = zeros(NumberWindows,1);
    for k = 1:NumberWindows
        NumberWindowsstore(k,1) = f(i,k) + b(i,k);
    end
    C(i,1) = 1/sum(NumberWindowsstore);
end
% Fourth and final step - Put all parameters together to solve for the
% intensities
Iclean = zeros(time,1);
for i = 1:time
    TempSum = zeros(NumberWindows,1);
    for k = 1:NumberWindows
        TempSum(k,1) = f(i,k)*C(i,1)*I_avg_f(i,k) + b(i,k)*C(i,1)*I_avg_b(i,k);
    end
    Iclean(i,1) = sum(TempSum);
end
IC = Iclean(2:time-1,1);
return
end




%{

time = length(rawTrace);



% start_index = max(windows) + M;
numberWindows = length(windows);
filteredTrace = zeros(1,length(rawTrace)-2*start_index);
% I_forward = zeros(numberWindows, length(rawTrace) - 2*start_index);
% I_reverse = zeros(numberWindows, length(rawTrace)-2*start_index);

% [numWindows, lengthTraj] = size(I_forward); 
                               % get the shape/size of I_forward
                               % used for iteration

                               
                           
I_avg_f = zeros(time, numberWindows);
I_avg_b = zeros(time, numberWindows);
                               
%}                               

%{

Pseudocoding the filter:
------------------------

iterate over the indexes of lengthTraj
    iterate over the index and values of windows
        calculate average of forward window
        calculate average of reverse window
    iterate over 

%}



%{
start_index = max(windows) + M
filtered_rawTrace = np.zeros(len(rawTrace)-2*start_index)
I_forward = np.zeros((len(windows), len(rawTrace) - 2*start_index))
I_reverse = np.zeros((len(windows), len(rawTrace) - 2*start_index))
for i in range(I_forward.shape[-1]):
    for j, k in enumerate(windows):
        I_forward[j,i] = forward_predictor(rawTrace, k, i+start_index-1)
        I_reverse[j,i] = reverse_predictor(rawTrace, k, i+start_index-1)
    for i in range(start_index, len(rawTrace) - start_index):
        f_k = np.zeros(len(windows))
        b_k = np.zeros(len(windows))
        for j, k in enumerate(windows):
            for m in range(M-1):
                f_k[j] += (rawTrace[i-m] - I_forward[j, i-m-start_index])**2
                b_k[j] += (rawTrace[i-m] - I_reverse[j, i-m-start_index])**2
            f_k[j] = f_k[j]**-p
            b_k[j] = b_k[j]**-p
        c = f_k.sum() + b_k.sum()

        for k in range(len(windows)):
            filtered_rawTrace[i-start_index] += f_k[k]*I_forward[k, i-start_index] + \
            b_k[k]*I_reverse[k, i-start_index]
        filtered_rawTrace[i-start_index] /= c
    return filtered_rawTrace
%}

   
  