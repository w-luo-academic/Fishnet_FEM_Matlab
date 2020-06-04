function [K] = stiff_Update_s(i_fail,K_old,connect,stiff_old,stiff_new)
%*************************************************************************
%
%           /====|    |=====|     |===|     |=====|    |=====|
%          |             |          |       |          |
%           \===\        |          |       |===|      |===|
%                |       |          |       |          |
%          |====/       |=|       |===|     |          |        _UPDATE
%
%*************************************************************************
%   update STIFFNESS MATRIX for m x n nacre system with SOFTENING
%                                                      ^^^^^^^^^^^
%       Input: i_fail = index of newly failed element
%              K_old = old stiffness matrix
%              connect = connectivity matrix (num._of_ele by 2)
%              stiff_old = old link stiffness to be replaced
%              stiff_new = new link stiffness to update
%
%       Output: K = New Global Stiffness Matrix
%
%       Warning: 1. m,n MUST be EVEN numbers!
%                2. This code works only with LINEAR ELASTIC material!
%
%   Wen Luo
%   8/22/2017
%*************************************************************************
%*************************************************************************

%   initialize global stiffness matrix K
    K = K_old;
    % node number 1 of Ke (element stiffness)
    n1 = connect(i_fail,1);
    % node number 2 of Ke
    n2 = connect(i_fail,2);
    % update Ke: k0 * [1,-1;-1,1]
    K(n1,n1) = K(n1,n1) - stiff_old + stiff_new;
    K(n1,n2) = K(n1,n2) + stiff_old - stiff_new;
    K(n2,n1) = K(n2,n1) + stiff_old - stiff_new;
    K(n2,n2) = K(n2,n2) - stiff_old + stiff_new;
end