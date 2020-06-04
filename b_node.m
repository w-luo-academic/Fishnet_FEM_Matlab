function [b_n]=b_node(m,n)
%*************************************************************************
%
%           ====\        ==   ==      ===     ====       =====|        
%           |    |       | \   |    /     \    |    \    |         
%           |===|        |  \  |   |       |   |     |   |==|      
%           |    |       |   \ |    \     /    |    /    |    
%           ====/  ----  ==   ==      ===     ====       =====|
%
%*************************************************************************
%   BOUNDARY NODE MATRIX for an m x n nacre system
%
%       Input: m = number of (zig-zag) rows (elements)
%              n = number of (zig-zag) columns (elements)
%
%       Output: b_n = boundary node matrix (m/2 by 2):
%               b_n(:,1) = LEFT boundary nodal numbers
%               b_n(:,2) = RIGHT boundary nodal numbers
%
%       Warning: m,n MUST be EVEN numbers!
%
%   Wen Luo
%   2/2/2017
%*************************************************************************
%*************************************************************************
    % total number of independent rows
    tot_n_row = m / 2;
    % initialize b_n matrix
    b_n = zeros(tot_n_row,2);
    % jump of nodal numbers between neiboring independent rows
    delta = n + 1;
    % evaluation loop
    for ii = 1:tot_n_row
       % left boundary
       b_n(ii,1) = 1 + (ii-1) * delta;
       % right boundary
       b_n(ii,2) = ii * delta;
    end
end