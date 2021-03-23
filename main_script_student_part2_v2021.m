clc;clear;
%LOADING DATASET
load Dataset_A.mat;

%HERE COME USER-DEFINED PARAMETERS OF GN algorithm (using default options is recommended)
eps_tol=10^-5;  %stopping criterion for max(abs(delta_x))
Max_iter=100; %maximum number of GN iterations
H_decoupled=0; %use full Jacobian matrix
H_sparse=1; %use matrix in a sparse format
linsolver=2;  %use Cholesky


%% DATA PREPARATION 
%getting Jacobian sparsity pattern (done only once)
[topo, ind_meas, N_meas]=f_meas_indices_and_H_sparsity_pattern_v2021(topo, meas);

%obtaining the admittance matrix in the compressed row form (done only once)
Y_bus = f_Y_bus_compressed_v2021(topo);



%constructing vector z, which is a vector of available measurements
%NOTE: different types of measurements should be in the same order as in the handout
%STUDENT CODE 1
z=???


%constructing the matrix of weights W and its square root Wsqrt using standard deviations provided in structure 'meas'
%NOTE: matrices W and Wsqrt must be constructed in sparse format!!!
%STUDENT CODE 2
W=???
Wsqrt=???




%% INITIAL SE SOLUTION
fprintf('INITIAL STATE ESTIMATION: \n')
%choosing the initial point
%STUDENT CODE 3
V=???
theta=???

%calling the GN algorithm and recording initial solution of SE (for
%subsequent plotting of voltage magnitudes and angles)
[ V_SE0, theta_SE0, ~, ~, convergence ] = f_SE_NR_algorithm_v2021 ( V, theta, topo, Y_bus, z, W, Wsqrt, ...
    ind_meas, N_meas, eps_tol, Max_iter, H_decoupled, H_sparse, linsolver );



%% BAD DATA DETECTION
%STUDENT CODE 4




%% BAD DATA CORRECTION
%STUDENT CODE 5





