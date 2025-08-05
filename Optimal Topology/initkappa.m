function [kappa_init] = initkappa(Omega_train, beta1, y_train)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  24/07/2025  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Developer: Sahand Tangerami %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function generates the initial \kappa for the optimization.    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:                                                              %
%       Omega_train:  State Matrix of The Reservoir Computer(Training)%
%       beta1:        Regularization parameter                        %
%       y_train:      Target Training Signal                          %
% Output:                                                             %
%       kappa_init:     Initial Guess of the Reservoir Output Weights %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    tol = 1e-2;
    kappa_init = (pinv(Omega_train' * Omega_train + beta1 * eye(size(Omega_train,2)), tol) * Omega_train' * y_train)';
end