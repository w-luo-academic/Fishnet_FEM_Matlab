function [v_new] = n_sort(v,b_n,i_n)
%*************************************************************************
%
%        =     =       /====|       ===       |===\     |=====| 
%        | \   |      |           /     \     |    |       |
%        |  \  |       \===\     |       |    |===/        |
%        |   \ |            |     \     /     |   \        |
%        =     = ____ |====/        ===       =    =      |=|
% 
%*************************************************************************
%   REARRANGE THE VECTOR IN SPLIT ORDER INTO ORIGINAL ORDER
%                        
%       Input: v = vector in split order (Dirichlet, Neumann)
%              b_n = boundary node number matrix
%              i_n = internal node number vector
%
%       Output: v_new = vector in original order (1,2,3,4,5,...)
%
%
%   Wen Luo
%   2/4/2017
%*************************************************************************
%*************************************************************************  
    % total number of nodes
    num_n = length(v);
    % initialize v_new
    v_new = zeros(num_n,1);
    % construct index vector
    ind_v = [b_n(:,1);b_n(:,2);i_n];
    % reorder the disp. vector
    for ii = 1: num_n
        index = ind_v(ii);
        v_new(index) = v(ii); 
    end
end