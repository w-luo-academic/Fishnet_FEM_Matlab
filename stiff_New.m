function [K] = stiff_New(m,n,connect,elas,e_area,e_len)
%*************************************************************************
%
%           /====|    |=====|     |===|     |=====|    |=====|
%          |             |          |       |          |
%           \===\        |          |       |===|      |===|
%                |       |          |       |          |
%          |====/       |=|       |===|     |          |          _ NEW
%
%*************************************************************************
%   STIFFNESS MATRIX for a linear elastic - brittle m x n nacre system
%
%       Input: m = number of (zig-zag) rows (elements)
%              n = number of (zig-zag) columns (elements)
%              connect = connectivity matrix (num._of_ele by 2)
%              elas = elastic modulus
%              e_area = element cross-section area
%              e_len = element length
%
%       Output: K = Global Stiffness Matrix
%
%       Warning: 1. m,n MUST be EVEN numbers!
%                2. This code works only with LINEAR ELASTIC material!
%                3. This stiffness matrix is for the INITIAL step ONLY!
%
%   Wen Luo
%   5/16/2017
%*************************************************************************
%*************************************************************************
%   total number of nodes
    tot_num_n = ( n + 1 ) * m / 2 + n / 2;
%   total number of elements
    tot_num_e = m * n;
%   modulus multiplier
    k0 = elas * e_area / e_len;
%   initialize global stiffness matrix K
    K = zeros(tot_num_n);
%   element loop
    for ii = 1:tot_num_e
        % node number 1 of Ke (element stiffness)
        n1 = connect(ii,1);
        % node number 2 of Ke
        n2 = connect(ii,2);
        % update Ke: k0 * [1,-1;-1,1]
        K(n1,n1) = K(n1,n1) + k0;
        K(n1,n2) = K(n1,n2) - k0;
        K(n2,n1) = K(n2,n1) - k0;
        K(n2,n2) = K(n2,n2) + k0;
    end
end