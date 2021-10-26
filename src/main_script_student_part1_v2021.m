clear;
%LOADING THE DATA
load GoodMeasurement_14_bus.mat;
%load GoodMeasurement_1354_bus.mat;

num_iter = 1; %Number of Iterations for the computational time average

%HERE COME USER-DEFINED PARAMETERS OF GN algorithm
eps_tol=10^-5;  %stopping criterion for max(abs(delta_x))
Max_iter=100; %maximum number of GN iterations
H_decoupled=0; %if 0, full H and 'normal' GN are used. If 1, decoupled H and fast decoupled GN are used
H_sparse=0; %if 1, H is created as a sparse matrix, if 0 - as a dense matrix
linsolver=4;  %1 - matrix inverse, 2 - Cholesky, 3 - QR, 4 - Hybrid

%SOLTUION ALGORITHM STARTS HERE
%getting Jacobian sparsity pattern (done only once)
[topo, ind_meas, N_meas]=f_meas_indices_and_H_sparsity_pattern_v2021(topo, meas);

%obtaining the admittance matrix in the compressed row form (done only once)
Y_bus = f_Y_bus_compressed_v2021(topo);

%constructing vector z, which is a vector of available measurements (hint: use structure ind_meas)
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
 
% %constructing the matrix of weights W and its square root Wsqrt using standard deviations provided in structure 'meas'
% %NOTE: matrices W and Wsqrt must be constructed in sparse format!!!
% %STUDENT CODE 2
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
 
% %choosing the initial point (use flat start)
% %STUDENT CODE 3
V = ones(topo.nBus,1);
theta= zeros(topo.nBus,1); % bus 1 is used as reference
 
% %TASK 1: SOLVING SE PROBLEM WITH GN ALGORITHM
% %STUDENT CODE 4
% %NOTE: use the following function (the description of its inputs and
% outputs can be found inside the function)
% [ V, theta, eps_all, time, convergence, it_num ] = f_SE_NR_algorithm_v2021 ( V, theta, topo, Y_bus, z, W, Wsqrt, ...
%             ind_meas, N_meas, eps_tol, Max_iter, H_decoupled, H_sparse, linsolver );         
tot_time = 0;
for n = 1 : num_iter
    V = ones(topo.nBus,1);
    theta= zeros(topo.nBus,1);
    [ V, theta, eps_all, time, convergence, it_num ] = f_SE_NR_algorithm_v2021 ( V, theta, topo, Y_bus, z, W, Wsqrt, ...
             ind_meas, N_meas, eps_tol, Max_iter, H_decoupled, H_sparse, linsolver );
    tot_time = tot_time + time;
end;
avg_time = tot_time/num_iter;
fprintf('Total time: %.6f s\n', tot_time);
fprintf('Average time: %.6f ms\n', avg_time*1000);
fprintf('Number iterations: %d \n', it_num);


% %TASK 2: COMPUTING CONDITION NUMBERS AND DENSITY FACTORS OF MATRICES
% %STUDENT CODE 5
% Set H_sparse = 1
%Density factor of H
[ H ] = f_measJac_H_v2021( V, theta, Y_bus, topo, ind_meas, N_meas, H_decoupled, H_sparse);
DF_H = nnz(H)/numel(H)*100;
%Density factor of G
G = H'*W*H;
DF_G = nnz(G)/numel(G)*100;
%Density factor of inverse of G
Ginv = inv(G);
DF_Ginv = nnz(Ginv)/numel(Ginv)*100;
% uncomment next 3 lines for the density factor
%Cholesky decompostion and density factor of L
% [L,p,S] = chol(G);
% DF_L = nnz(L)/numel(L)*100;
%QR decomposition and densitiy factor of R an Q
[Q,R,e] = qr(H,0);
DF_R = nnz(R)/numel(R)*100;
DF_Q = nnz(Q)/numel(Q)*100;

%Set H_sparse = 0
%Condition number of G
cn_G = cond(G);
%Condition number of H & R
cn_R = cond(R);
 
% %TASK 3: COMPARING PERFORMANCE OF FULL/DECOUPLED GN
% %STUDENT CODE 6
% %NOTE: use the following function (the description of its inputs and
% %outputs can be found inside the function)
H_decoupled = 1;

tot_time = 0;
for n = 1 : num_iter
    V = ones(topo.nBus,1);
    theta= zeros(topo.nBus,1);
    [ V, theta, eps_all, time, convergence, it_num] = f_SE_NR_algorithm_v2021 ( V, theta, topo, Y_bus, z, W, Wsqrt, ...
            ind_meas, N_meas, eps_tol, Max_iter, H_decoupled, H_sparse, linsolver );
    tot_time = tot_time + time;
end;
avg_time = tot_time/num_iter;
fprintf('Total time: %.6f s\n', tot_time);
fprintf('Average time: %.6f ms\n', avg_time*1000);
fprintf('Number iterations: %d \n', it_num);