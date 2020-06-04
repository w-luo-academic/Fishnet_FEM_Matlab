function [sig_0] = r_sig_G(num_e,mean,std_dev)
%*************************************************************************
%
%         |===\          /====|    |=|     |===|          |===|
%         |    |        |           |     /              /
%         |===/          \===\      |    |     ==|      |     ==|
%         |   \               |     |     \     /        \     /
%         =    =  ____  |====/     |=|     |===|   ____   |===|
%  
%*************************************************************************
%   RANDOM STRENGTH GENERATOR (GAUSSIAN)
%                        
%       Input: num_e = total number of elements
%              mean = mean of a single element
%              std_dev = standard deviation of a single element
%
%       Output: sig_0 = vector of random strength
%
%
%   Wen Luo
%   2/4/2017
%*************************************************************************
%*************************************************************************  
rng('shuffle');
sig_0 = std_dev.*randn(num_e,1) + mean;
end