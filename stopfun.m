%Output function to stop solver if it takes too long

function stop = stopfun(x,optimValues,state)
T=5; %Seconds
stop = toc>T;