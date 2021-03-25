clc;clear;
%LOADING DATASET
load Dataset_B.mat;

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
fields = fieldnames(meas);
z=[];
for i = 1:numel(fields)
    if ~contains(fields{i},'warning') && ~contains(fields{i},'std')
        z = [z;meas.(fields{i})];
        z = z(~isnan(z));
    end;
end;

%constructing the matrix of weights W and its square root Wsqrt using standard deviations provided in structure 'meas'
%NOTE: matrices W and Wsqrt must be constructed in sparse format!!!
%STUDENT CODE 2
fields = fieldnames(meas);
all_std=[];
all_msr=[];

% Get all std deviation elements + measurements in two different vectors
for i = 1:numel(fields)
  if ~contains(fields{i},'warning')
      if contains(fields{i},'std')
        all_std = [all_std; meas.(fields{i})];
      else
        all_msr = [all_msr; meas.(fields{i})];
      end;
  end;
end

% Fill all 'NaN' elements (from measurement vector) in the std deviation
% vector
for i = 1 : size(all_msr,1)
    if isnan(all_msr(i))
        all_std(i) = NaN;
    end;
end;

% Remove all nodes which are not measured (filter NaN)
msr_std = all_std(~isnan(all_std));

% Diagonal elements of W is given by 1/(std^2)
W = sparse(inv(diag(msr_std.^2)));
Wsqrt = sqrt(W);


%% INITIAL SE SOLUTION
fprintf('INITIAL STATE ESTIMATION: \n')
%choosing the initial point
%STUDENT CODE 3
V = ones(topo.nBus,1);
theta= zeros(topo.nBus,1); % bus 1 is used as reference

%calling the GN algorithm and recording initial solution of SE (for
%subsequent plotting of voltage magnitudes and angles)
[ V_SE0, theta_SE0, ~, ~, convergence ] = f_SE_NR_algorithm_v2021 ( V, theta, topo, Y_bus, z, W, Wsqrt, ...
    ind_meas, N_meas, eps_tol, Max_iter, H_decoupled, H_sparse, linsolver );



%% BAD DATA DETECTION
%STUDENT CODE 4
[ H ] = f_measJac_H_v2021(V_SE0, theta_SE0, Y_bus, topo, ind_meas, N_meas, H_decoupled, H_sparse);
[ h ] = f_measFunc_h_v2021( V_SE0, theta_SE0, Y_bus, topo, ind_meas, N_meas);

r = z - h;

G = H' * W * H;  
K = H * (G \ H') * W;
S = round(diag(eye(size(K,1)) - K),12); %remove numerical instabilities which lead to negative values
cov_R = S./(diag(W));

r_norm = abs(r)./sqrt(cov_R); %there are infinity values because some zeroes are created in the covariance matrix(?)



%% BAD DATA CORRECTION
%STUDENT CODE 5





