
function IC = NoiseReduc(NoisyTrace)
%
%Coded by Nigel Reuel on 8.31.2010
%
%Adaptation of Chung and Kennedy (1991) forward-backward non-linear
%filtering technique to reduce the noise of our
%fluorescent signal in an attempt to resolve the single molecule events
%amidst the noise.
%
%Input: Takes supplied noisy trace (NoisyTrace)
%Output: Returns cleaned trace without the first and last data point
%
time = length(NoisyTrace);
%
%Set algorithm parameters (discussion on this in Bioinformatics paper - Supplement A):
%
% ^^^ NOTE: User can change following parameter to improve fit ^^^
%
window_list = [4 8 16 32]; %Sampling window length
number_of_windows = length(window_list); %Number of sampling windows
M = 10; %Weighting parameter length
P = 40; %Parameter that effects sharpness of steps

% ——- The Forward-Backward non-linear Algorithm ———-

%First - Run each time step and store the average forward and backward
%predictors for each of the sampling windows

I_avg_f = zeros(time,number_of_windows);
I_avg_b = zeros(time,number_of_windows);
for time_point = 1:time
%Run for each sampling window
    for window = 1:number_of_windows
    %Solve for the average forward predictor (w/ some
    %logic to help solve the beginning of the trace
    window_size = window_list(window);
        if time_point == 1
            I_avg_f(time_point,window_size) = NoisyTrace(1,1);
        elseif time_point - window_size - 1 < 0
            I_avg_f(time_point,window_size) = ...
                sum(NoisyTrace(1:time_point-1,1))/time_point;
        else
            end_point = time_point - window_size;
            start_point = time_point - 1;
            I_avg_f(time_point,window_size) = ...
                sum(NoisyTrace(end_point:start_point,1))/window_size;
        end
%Now do the same for the backward predictor
        if time_point == time
            I_avg_b(time_point,window_size) = NoisyTrace(time_point,1);
        elseif time_point + window_size > time
            sw = time - time_point;
            I_avg_b(time_point,window_size) = ...
                sum(NoisyTrace(time_point+1:time,1))/sw;
        else
            end_point = time_point + window_size;
            start_point = time_point + 1;
            I_avg_b(time_point,window_size) = ...
                sum(NoisyTrace(start_point:end_point,1))/window_size;
        end
    end
end

%Second - Solve for non-normalized forward and backward weights:
f = zeros(time,number_of_windows);
b = zeros(time,number_of_windows);
for time_point = 1:time
    for window_size = 1:number_of_windows
        Mstore_f = zeros(M,1);
        Mstore_b = zeros(M,1);
        for j = 0:M-1
            t_f = time_point - j;
            t_b = time_point + j;
            if t_f < 1
                Mstore_f(j+1,1) = ...
                    (NoisyTrace(time_point,1) - I_avg_f(time_point,window_size))^2;
            else
                Mstore_f(j+1,1) = ...
                    (NoisyTrace(t_f,1) - I_avg_f(t_f,window_size))^2;
            end
            if t_b > time
                Mstore_b(j+1,1) = ...
                    (NoisyTrace(time_point,1) - I_avg_b(time_point,window_size))^2;
            else
                Mstore_b(j+1,1) = ...
                    (NoisyTrace(t_b,1) - I_avg_b(t_b,window_size))^2;
            end
        end
        f(time_point,window_size) = sum(Mstore_f)^(-P);
        b(time_point,window_size) = sum(Mstore_b)^(-P);
    end
end
% Third - Solve for vector of normalization factors for the weights
C = zeros(time,1);
for time_point = 1:time
    Kstore = zeros(number_of_windows,1);
    for window = 1:number_of_windows
        Kstore(window,1) = f(time_point,window) + b(time_point,window);
    end
    C(time_point,1) = 1/sum(Kstore);
end
% Fourth and final step - Put all parameters together to solve for the
% intensities
Iclean = zeros(time,1);
for time_point = 1:time
    TempSum = zeros(number_of_windows,1);
    for window = 1:number_of_windows
        TempSum(window,1) = f(time_point,window)*C(time_point,1)*I_avg_f(time_point,window) + b(time_point,window)*C(time_point,1)*I_avg_b(time_point,window);
    end
    Iclean(time_point,1) = sum(TempSum);
end
IC = Iclean(2:time-1,1);
return
end
