%ALGORITHM Put the algorithm you want to use here
%   Change the algorithm below to test your desired algorithm
%   data.x = input conditions
%   data.y = responses
%   lowerbounds = lower bounds of conditions
%   upperbounds = upper bounds of conditions

    Opt = TSEMO_options();  
    conditions = TSEMO(data.x, data.y, lowerbounds, upperbounds, Opt); 