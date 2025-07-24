function [y_train,y_test] = signal_creation(omega,t, b, phi, len_trans, len_train)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  24/07/2025  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Developer: Sahand Tangerami %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Function Generates The Training and Test Signals.              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:                                                              %
%       omega:      Input Signal Frequencies                          %
%       t:          Total Simulation Time Vector                      %
%       b:          Amplitudes of Output Signals                      %
%       phi:        Output Signal Phase Change                        %
%       len_trans:  Training Signal Length                            %
%       len_train:  Test Signal Length                                %
% Output:                                                             %
%       y_train:         Training Signal                              %
%       y_test:         Training Signal                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    y = generate_output(t, b, omega, phi);
    
    y_train = y(len_trans+1:len_trans+len_train,:);
    y_test = y(len_trans+len_train + 1:end,:);
end