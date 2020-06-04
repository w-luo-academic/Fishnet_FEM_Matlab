function [conn_R] = conn_row(start, n)
%*************************************************************%
%   Connectivity Matrix for a single INDEPENDENT row
%
%       Input: start = initial node number
%              n     = total number of elements of that row
%
%       Output: conn_R = connectivity matrix
%
%   Wen Luo
%   1/20/2017
%
%*************************************************************%
    conn_R = zeros(n,2);
    for ii = 1:n
        conn_R(ii,:) = [start+ii-1,start+ii];
    end
end