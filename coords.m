function [coord] = coords(m,n)
%*************************************************************************
%
%                              coords_2D
%                                                        
%*************************************************************************
%   COORDINATES MATRIX for an m x n nacre system (2D)
%
%       Input: m = number of (zig-zag) rows (elements)
%              n = number of (zig-zag) columns (elements)
%
%       Output: coord = coordinates matrix of nacre system:
%               coord(i,j) = the x_j coord. of node i
%
%       Warning: m,n MUST be EVEN numbers!
%
%   Wen Luo
%   7/10/2017
%   -----------------------------------------------------------
%   ALGORITHM:
%                      Node Number         
%                   2   4   6   8   10     
%   i=1            / \ / \ / \ / \ / \     
%                 1   3   5   7   9   11   
%      
%                   13  15  17  19  21    
%   i=2            / \ / \ / \ / \ / \    
%                 12  14  16  18  20  22
%   ...
%
%   i=m/2          / \ / \ / \ / \ / \
%
%   last row       \ / \ / \ / \ / \ /
%          
%   Each bar length is 1x10^-3 mm (1 um) and Node 1 has coord x = 0.
%*************************************************************************
%*************************************************************************
%   total number of nodes:   last row has n/2 nodes
    num_node = (n+1) * m / 2 + n / 2;
%   initialize matrix coord
    coord = zeros (num_node,2);
    
%   assign nodal coordinates:

%-------------------------------------------------------------------------
%               First m/2 Independent Rows
%-------------------------------------------------------------------------
%   independent row loop
    for ii = 1:m/2
        % total number of nodes on each row
        delta = n + 1;
        % node # of the first one at current row
        start_node = (ii-1) * delta + 1;
        % node # of the last one at current row
        end_node = ii * delta;
        % beginning y coord
        y_start = (2*ii-2) * 0.5;
        % end y coord
        y_end = y_start + 0.5;
        % assign nodal coords for current row : from 0 to n
        coord(start_node:end_node,1) = 0:n;
        %
        jj = 1;
        while jj <= n/2
            coord(start_node+(2*jj-1),2) = y_start;
            coord(start_node+(2*jj-2),2) = y_end;
            coord(end_node,2) = y_end;
            jj = jj + 1;
        end
    end
%-------------------------------------------------------------------------
%               Last Row
%-------------------------------------------------------------------------
%   initialize last row coords vector
    coords_end = zeros(n/2,1);
%   evaluate last row coords
    for ii = 1:n/2
        coords_end(ii) = 2 * ii - 1;
    end
%   node # of first one of last row
    start_node = (n+1) * m / 2 + 1;
%   node # of last one of last row
    end_node = start_node + n / 2 - 1;
%   assign nodal coords for last row : n/2 odd numbers
    coord(start_node:end_node,1) = coords_end;
    coord(start_node:end_node,2) = m * 0.5;
end