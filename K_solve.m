function [u_2,f_R] = K_solve(K11,K12,K21,K22,u_1)
%*************************************************************************
%
%    =    =         /====|       ===       |        =      =     =====|
%    |   /         |           /     \     |         \    /      |
%    |==\           \===\     |       |    |          \  /       |==|
%    |   \               |     \     /     |           \/        |
%    =    =  ____  |====/        ===       |=====|     ==        =====|
% 
%*************************************************************************
%   REARRANGE THE VECTOR IN SPLIT ORDER INTO ORIGINAL ORDER
%                        
%       Input: K11,K12,K21,K22 = block stiff matrices
%              u_1 = Dirichlet B.C. nodal displacement vector (boundary)
%
%       Output: u_2 = Neumann B.C. nodal displacement vector (internal)
%               f_R = Reaction forces on Dirichlet B.C. nodes (f_1) (boundary)
%
%
%   Wen Luo
%   2/4/2017
%*************************************************************************
%*************************************************************************  
    % check if K22 is singular
    if det(K22) <= 0
        %error('K22 is SINGULAR, structure FAILED COMPLETELY!')
    end
    %----------------------------------------------------
    %     LATER ADD STOP VARIABLES OR STOP CONDITIONS
    %----------------------------------------------------
    
    % det(K22) is still positive definite:
    u_2 = - K22 \ (K21 * u_1);
    f_R = [K11, K12] * [u_1; u_2];
end