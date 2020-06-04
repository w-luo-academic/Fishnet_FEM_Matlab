function [stress,f_1,u] = kernel_New(K,k_vector,connect,area,length,b_n,i_n,n_1,n_2,num_n,num_e,u_1)
%*************************************************************************
%
%        =    =    =====|    |===\     =     =     =====|    =                   
%        |   /     |         |    |    | \   |     |         |     
%        |==\      |==|      |===/     |  \  |     |==|      |  
%        |   \     |         |   \     |   \ |     |         |             
%        =    =    =====|    =    =    =     =     =====|    ====== _ NEW                          
% 
%*************************************************************************
%   FEM Kernel of MXN NACRE TRUSS ELASTIC SYSTEM
%                        
%       Input: K = stiffness matrix
%              k_vector = stiffness status flag vector (num._of_ele by 1)
%              connect = connectivity matrix (num._of_ele by 2)
%              area = cross section area (mm^2)
%              length = length of element (mm)
%              b_n = boundary node number matrix
%              i_n = internal node number vector
%              n_1 = number of Dirichlet B.C nodes
%              n_2 = number of Neuman B.C nodes
%              num_n = total number of nodes
%
%      ***>>> Loading:
%              u_1 = Dirichlet B.C. nodal displacement vector (boundary)
%
%       Output: stress = stress vector of each element
%               u = nodal displacement vector (original order 1,2,3,...)
%               f_1 = Reaction forces on Dirichlet B.C. nodes (boundary)
%
%
%   Wen Luo
%   2/4/2017
%
%   8/22/2017 updated: F_flag replaced by softening status k_vector
%
%*************************************************************************
% COMMENT on indexing of u and f:
% Initially nodal numbers are assign based on the convention given in
% finction "conn". In order to solve for unknown displacements and reaction
% forces, the stiffness matrix K is rearranged into 4 block matrices, so
% the nodal arrangement is changed from original one to:
%
%           [left boundary; right boundary; internal]
%
% After all unknowns are solved, the nodal displacement u and force vector
% f is in the new order of arrangement. So function "n_sort" is used to
% convert them back to the original order:
%                       [1,2,3,4,5,...]
%*************************************************************************  

    % split K into block matrices for solving the system
    [K11,K12,K21,K22] = split(K,b_n,i_n,n_1,n_2,num_n);

    % solve for u_2 and f_R
    [u_2,f_1] = K_solve(K11,K12,K21,K22,u_1);

    % construct disp. vector in SPLIT order
    u_s = [u_1;u_2];

    % rearrange u to the original order
    u = n_sort(u_s,b_n,i_n);

    % stress and strain vector
    strain = zeros(num_e,1);
    stress = strain;
    for ii = 1:num_e
        ind_1 = connect(ii,1);
        ind_2 = connect(ii,2);
        strain(ii) = ( u(ind_2) - u(ind_1) ) / length;
        stress(ii) = strain(ii) * k_vector(ii) * length / area;
    end
end