function [K11,K12,K21,K22]=split(K,b_n,i_n,n_1,n_2,n_node)
%*************************************************************************
%
%           /====|    |====\     |          |===|     |=====|    
%          |          |     |    |            |          |           
%           \===\     |====/     |            |          |      
%                |    |          |            |          |          
%          |====/     =          |=====|    |===|       |=|          
%
%*************************************************************************
%   SPLIT THE STIFFNESS MATRIX K INTO 4 SEPERATE BLOCK MATRICES
%                        
%                         / K11 | K12 \ / u_0 \
%                   K u = |  -------  | | --- |
%                         \ K21 | K22 / \ u_1 /
%
%           u_0: nodal disp. with Dirichlet B.C.
%           u_1: nodal disp. with Neumann B.C
%
% ------------------------------------------------------------------------
%       Input: K = Global Stiffness Matrix
%              b_n = boundary node number matrix
%              i_n = internal node number vector
%              n_1 = number of Dirichlet B.C nodes
%              n_2 = number of Neuman B.C nodes
%              n_node = total number of nodes
%
%       Output: [K11,K12,K21,K22]
%
%       Warning: 1. m,n MUST be EVEN numbers!
%                2. This code works only with LINEAR ELASTIC material!
%
%   Wen Luo
%   2/3/2017
%*************************************************************************
%*************************************************************************    
    %------ K11 ------
    % initialize K11 and K_temp
    K11 = zeros(n_1);
    K_temp = zeros(n_1,n_node);
    % pick rows
    for ii = 1:n_1
        index = b_n(ii);
        K_temp(ii,:) = K(index,:);
    end
    % pick colomns
    for ii = 1:n_1
        index = b_n(ii);
        K11(:,ii) = K_temp(:,index);
    end
    
    %%------ K12 ------
    % initialize
    K12 = zeros(n_1,n_2);
    % pick colomns (rows are the same with K11_temp)
    for ii = 1:n_2
        index = i_n(ii);
        K12(:,ii) = K_temp(:,index);
    end
    
    %%------ K22 ------
    % initialize K22 and K_temp
    K22= zeros(n_2);
    K_temp = zeros(n_2,n_node);
    % pick rows
    for ii = 1:n_2
        index = i_n(ii);
        K_temp(ii,:) = K(index,:);
    end
    % pick colomns
    for ii = 1:n_2
        index = i_n(ii);
        K22(:,ii) = K_temp(:,index);
    end
    
    %%------ K21 ------
    % transpose K12 to get K21
    K21 = K12';
end