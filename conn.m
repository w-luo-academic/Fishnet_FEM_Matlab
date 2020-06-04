function [connect] = conn (m,n)
%*************************************************************************
%
%                  ===      ===      ==   ==    ==   ==  
%                ||       ||   ||    ||\  ||    ||\  ||
%                ||       ||   ||    || \ ||    || \ ||
%                ||       ||   ||    ||  \||    ||  \||
%                  ===      ===      ==   ==    ==   ==
%
%*************************************************************************
%   CONNECTIVITY MATRIX for an m x n nacre system
%
%       Input: m = number of (zig-zag) rows (elements)
%              n = number of (zig-zag) columns (elements)
%
%       Output: connect = connectivity matrix of nacre system:
%               connect(i,j) = the j'th node number of element i
%
%       Warning: m,n MUST be EVEN numbers!
%
%   Wen Luo
%   1/20/2017
%   -----------------------------------------------------------
%   LABELING CONVENTION:
%
%         Node Number         |         Element Number
%      2   4   6   8   10     |
%     / \ / \ / \ / \ / \     |   1, 2, 3, 4, 5, 6, 7, 8, 9,10
%    1   3   5   7   9   11   |
%     \ / \ / \ / \ / \ /     |  11,12,13,14,15,16,17,18,19,20
%      13  15  17  19  21     |
%     / \ / \ / \ / \ / \     |  21,22,23,24,25,26,27,28,29,30 
%    12  14  16  18  20  22   |    (Last row will be deleted)
%
%   ------------------------------------------------------------
%   ALGORITHM:
%                   n colomns
%          ----------------------------
%          | /\/\/\/\/\/\/\/\/\/\/\/\ |    Independent
%          | \/\/\/\/\/\/\/\/\/\/\/\/ |    Dependent
%   m rows | /\/\/\/\/\/\/\/\/\/\/\/\ |    Independent
%          | \/\/\/\/\/\/\/\/\/\/\/\/ |    Dependent
%          |           ...            |    ...
%          | /\/\/\/\/\/\/\/\/\/\/\/\ |    Independent (To be deleted)
%          ----------------------------
%
%   Entries of connect matrix for independent rows are the same (similar) 
%   as the 1-D case; Those for dependent rows are formed by the two 
%   neighboring independent rows: 
%
%           A11  A12
%   n rows  A21  A22
%           A31  A32        
%           ...              
%          ----------      where A stands for corresponding
%           A11  B12       entries from the same location of
%   n rows  B21  A22       the PREVIOUS independent block; 
%           A31  B32       B stands for those from the LATTER 
%           ...            independent block.
%          ----------
%           B11  B12
%   n rows  B21  B22
%           B31  B32   
%           ...
%     
%*************************************************************************
%*************************************************************************
    connect = zeros((m+1)*n,2);
    % element number of first element of the last row
    last_1 = m*n+1;
    % last row nodal number
    lastR_node = linspace(m*n/2+m/2+1, m*n/2+m/2+n/2, n/2);
    
    % ------------------------------------------------
    %               INDEPENDENT ROWS: 
    % ------------------------------------------------
    for ii = 1:m/2+1 % one more extra block!
        % number of nodes in each row: n+1
        % beginning nodal number: (i-1)(n+1)+1
        start_node_N = (ii-1) * (n+1) + 1;
        
        % >>>
        % 2n(i-1) is the No. of ele's in front of the beginning of i'th
        % independent row.
        % beginning element number (beginning row index of connect):
        start_ele_N = 2 * n * (ii-1) + 1;
        % end element number (end row index of connect)
        end_ele_N = 2 * n * (ii-1) + n;
        connect(start_ele_N:end_ele_N,:) = conn_row(start_node_N,n);
    end
    % >>>
    % The last block is used to construct the second last block, so only 
    % half of its entries are usful (1 entry useful for each row and every 
    % 2 rows have the same useful entry), so it suffices to evaluate each 
    % 4x4 block by the constant matrix: lastR_node(ii)*[1,1;1,1]
    %
    % modify the last block 
    for ii = 1:n/2
        ind = last_1+2*(ii-1);
        %>>>
        % ind:   index of odd rows of the last block
        % ind+1: index of even rows of the last block
        connect(ind:ind+1,:) = lastR_node(ii)*[1,1;1,1];
    end
    
    % ------------------------------------------------
    %                 DEPENDENT ROWS: 
    % ------------------------------------------------
    for ii = 1:m/2
        % >>>
        % both start and end element numbers of dependent rows have
        % difference n with their independent counterparts.
        start_ele_N =  2* n * (ii-1) + n + 1;
        end_ele_N = 2 * n * (ii-1) + 2 * n;
        % evaluate connectivity matrix
        for jj = 1:n/2
            ind = 2 * jj - 1;
            % odd rows
            connect(start_ele_N+ind-1,1) = connect(start_ele_N+ind-1-n,1);
            connect(start_ele_N+ind-1,2) = connect(start_ele_N+ind-1+n,2);
            % even rows
            connect(start_ele_N+ind,1) = connect(start_ele_N+ind+n,1);
            connect(start_ele_N+ind,2) = connect(start_ele_N+ind-n,2);
        end
    end
    % delete the last block after the second last block is evaluated
    connect = connect(1:m*n,:);
end