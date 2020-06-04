function [inner_n] = inner_node(b_n,connect,m,n)
%*************************************************************************
%
%             O          ==   ==      ===     ====       =====|        
%             |          | \   |    /     \    |    \    |         
%             |          |  \  |   |       |   |     |   |==|      
%             |          |   \ |    \     /    |    /    |    
%            ===  -----  ==   ==      ===     ====       =====|
%
%*************************************************************************
%   INTERNAL NODE MATRIX for an m x n nacre system
%
%       Input: b_n = boundary node matrix (m/2 by 2)
%              connect = conectivity matrix
%              m = number of (zig-zag) rows (elements)
%              n = number of (zig-zag) columns (elements)
%
%       Output: inner_n = internal node vector: list of inner node numbers
%               
%
%       Warning: m,n MUST be EVEN numbers!
%
%   Wen Luo
%   2/3/2017
%*************************************************************************
%*************************************************************************
    % total num. of nodes
    num_n = (n+1) * m/2 + n/2;
    % number of prescribed disp. nodes (b_node)
    n_1 = numel(b_n);
    % number of internal nodes
    n_2 = num_n - n_1;
    % initialization
    inner_n = zeros(n_2,1);
    % independent rows loop
    for ii = 1:m/2
       % start index in inner_n vector
       start_ind = 1 + (ii-1) * (n-1);
       % end index in inner_n vector
       end_ind = ii * (n-1);
       % start node number of that independent row
       start_n = 2 + (ii-1) * (n+1);
       % end node number of that independent row
       end_n = n + (ii-1) * (n+1);
       % evaluate inner_n for that independent row
       inner_n(start_ind:end_ind) = start_n:end_n;
    end
    % last row
    % node number of the last one in internal nodes
    end_n = connect(end,1);
    % node number of the fisrt one of last row
    start_n = end_n - n/2 + 1;
    % evaluate inner_n for the last row
    inner_n(end - n/2 + 1:end) = start_n:end_n;
end