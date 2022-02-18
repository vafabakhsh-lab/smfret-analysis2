function  fwdAverage = forward_predictor(traj, N, i)

%{

This is assistant function for the chung-kennedy nonlinear filter. 
Source: https://doi.org/10.1016/0165-0270(91)90118-J

It calculates the average intensity in the forward window.

    Parameters
    ----------
    traj: numpy array
        trajectory
    N: int
        size of window
    i: int
        index
    Returns
    -------
    average intensity of forward window

%}
    
fwdAverage = sum(traj((-N:i-1))/N);


