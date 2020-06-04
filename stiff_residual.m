function [k_r] = stiff_residual(k0,sig_ori,sig_r,k_t)
%*************************************************************************
%
%                           STIFF_RESIDUAL
%
%*************************************************************************
%   Calculate the residual stiffness of a link in SOFTENING range
%                                                ^^^^^^^^^^^^^^^^^
%       Input: k0      = initial (original) link stiffness
%              sig_ori = original strength vector 
%              sig_r   = residual strength vector 
%              k_t     = tangential softening stiffness (NEGATIVE!!!)
%
%       Output: k_r = residual stiffness of the latest damaged link
%
%       Warning: 1. m,n MUST be EVEN numbers!
%                2. k_t is NEGATIVE!
%
%   Wen Luo
%   8/22/2017
%*************************************************************************
%*************************************************************************
    k_r = sig_r / (sig_ori/k0 - sig_ori/k_t + sig_r/k_t);
end