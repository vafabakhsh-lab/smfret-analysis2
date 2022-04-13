%{

Steps:
1. read in trace data
2. calculate the predictors for a given window size
    a. forward predictor
    b. backward predictor
3. calculate weights for each predictor
???
4. return filtered data


%}

% declare variables
window = 5;
M = 10;
p = 10;



% 'read' in trace data
NoisyTrace = [1 2 1 2 1 3 1 2 1 3 2 1 3 2 1 7 8 9 7 6 ...
    7 8 9 8 7 6 7 8 9 7 8 1 2 1 2 1 3 1 2 1 3 2 1 3 2 1];

NoisyTrace = NoisyTrace';


% determine starting point and end point of the trajectory because you
% can't calculate the average for i - window or i + window when you're at
% the end of the trajectory

starting_time = 1 + window;
ending_time = length(NoisyTrace) - (window + 1);


% calculate predictors

for current_time = starting_time:ending_time
    
    fwd_predictor(current_time) = ...
        mean(NoisyTrace(current_time - window : current_time - 1)) ;
    
    rev_predictor(current_time) = ...
       mean(NoisyTrace(current_time + 1 : current_time + window)) ;
        

end











