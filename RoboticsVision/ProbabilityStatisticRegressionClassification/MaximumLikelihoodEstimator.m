H=10; N=25;
mu_list=[0:0.01:1];
[max_L,max_mu]=MaxLikelihood(mu_list,H,N);

function [max_L,max_mu]=MaxLikelihood(mu_list,H,N)
    % compute the maximum likelihood value and its parameters
    %
    % Input:
    % - mu_list: a vector of possible parameter values, e.g. mu_list=[0:0.01:1]
    % - H: number of Heads
    % - N: number of total coin flips
    % Output:
    % - max_L: maximum likelihood value over all the parameters mu_list
    % - max_mu: a parameter value corresponding to the maximum likelihood value max_L
    
    % Compute a vector of likelihoods L_list
    % Every element in L_list should store a likelihood value
    % associated with its respective parameter from mu_list.
    L_list= mu_list.^H.*(1-mu_list).^(N-H);
    
    % Compute the maximum likelihood value by selecting a maximum
    % value from L_list
    [max_L, index] = max(L_list);
    
    % Find the parameter mu that corresponds to the maximum 
    % likelihood value
    max_mu= mu_list(index);  
end