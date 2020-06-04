function [sig_0] = r_sig_PG(num_e)
%*************************************************************************
%
%    |===\          /====|    |=|     |===|           |====\     |===|  
%    |    |        |           |     /                |     |   /       
%    |===/          \===\      |    |     ==|         |====/   |     ==|
%    |   \               |     |     \     /          |         \     / 
%    =    =  ____  |====/     |=|     |===|   ____    =          |===|
%  
%*************************************************************************
%   RANDOM STRENGTH GENERATOR (POWER LAW - GAUSSIAN)
%                        
%       Input: num_e = total number of elements
%
%       Output: sig_0 = vector of random strength
%
%
%   Wen Luo
%   2/6/2017
%*************************************************************************
%*************************************************************************  
rng('shuffle');
p0 = rand(num_e,1);
sig_0 = zeros(num_e,1);
a1 = 11.3379;
b1 = 38;
a2 = 8.83883;
b2 = 0.88883;
c2 = 0.00789331;
d2 = 0.503947;
for ii = 1:num_e
   if p0(ii) <= 0.0150364
       sig_0(ii) = 10*(p0(ii)/a1)^(1/b1);
   elseif p0(ii) > 0.0150364
       sig_0(ii) = (a2-erfinv( 1 - (p0(ii)+c2) / d2 )) / b2;
   end
end
end