function [ V, theta, eps_all, time, convergence, it_num ] = f_SE_NR_algorithm_v2021 ( V, theta, topo, Y_bus, z, W, Wsqrt, ...
    ind_meas, N_meas, eps_tol, Max_iter, H_decoupled, H_sparse, linsolver )
%**************************************************************************
%DESCRIPTION OF FUNCTION (VERSION 2021)
%**************************************************************************
%This function solves the State Estimation problem cast as a WLS problem
%using the Gauss-Newton (GN) method. The possible solution options include
%choosing different algorithms for solving a linear system at each
%iteration as well as using either a full GN method or fast decoupled GN
%method.

%INPUTS
%V           - vector of voltage magnitudes of ALL buses in the system
%theta       - vector of voltage phase angles for ALL buses in the system
%topo        - structure containing the topology of the system and
%              parameters of system elements. Used for computing H and h
%Y_bus       - structure that contains real and imag parts of admittances
%              of structurally nonzero elements in the admittance matrix
%z           - vector containing all available measurements in the system
%W           - diagonal matrix containing weights assigned to measurements
%Wsqrt       - diagonal matrix, elements are square roots of elements of W   
%ind_meas    - structure that for each type of measurement contains info
%              on what are the indexes of available measurements in a list
%              of all possible measurements for this type
%N_meas      - structure that contains the number of available measurements
%              for each type of measurements 
%eps_tol     - stopping criterion. GN should stop if max(|delta_x|)<eps_tol
%Max_iter    - stopping criterion, maximum number of GN iterations
%H_decoupled - control variable. If 0, full Jacobian is constructed and 
%              'full' GN algorithm is used. If 1, decoupled Jacobian is 
%              constructed and fast decoupled algorithm is used
%H_sparse    - control variable. If 1, matrix H will be in a sparse form. 
%              If 0, matrix H will be in a dense form
%linsolver   - control variable for choosing a linear system solver:
%              1 - matrix inverse of G
%              2 - Cholesky factorization of G (done by backslash)
%              3 - QR decomposition of H with minimum number of fill-ins
%              4 - Hybrid method (decompose H by QR, but only compute R)

%OUTPUTS
%V           - vector of voltage magnitudes of ALL buses in the system
%theta       - vector of voltage phase angles of ALL buses in the system
%eps_all     - vector, each element of which is equal to max(|delta_x|) at
%              a certain iteration. Shows the convergence pattern
%time        - total computation time for all iterations for solving linear
%              system. Time for computing H and h should NOT be included here
%convergence - if 1, GN has converged to a solution, otherwise 0

%**************************************************************************
%FUNCTION CODE
%**************************************************************************


%STUDENT CODE
%for computing matrix H and vector h, use the following functions:
eps_all = 0;
time = 0;
convergence = 0;
it_num = 0;

switch H_decoupled
    case 0
        for i = 1 : Max_iter
            [ H ] = f_measJac_H_v2021( V, theta, Y_bus, topo, ind_meas, N_meas, H_decoupled, H_sparse);
            [ h ] = f_measFunc_h_v2021( V, theta, Y_bus, topo, ind_meas, N_meas);


            tic;
            switch linsolver
                case 1 % Direct inverse
                    % 4. Compute G(x(k))
                    G = H' * W * H;  
                    % 5. Compute right hand side g(x(k))
                    g = -(H')*W*(z-h);
                    dx = -inv(G)*g;
                case 2 % Cholesky
                    % 4. Compute G(x(k))
                    G = H' * W * H;  
                    % 5. Compute right hand side g(x(k))
                    g = -(H')*W*(z-h);
                    dx = -G\g;
                case 3 % Orthogonal Decomposition
                    Htilda = Wsqrt*H;
                    [Q,R,e]=qr(Htilda,0);
                    b = Wsqrt*(z-h);
                    dx(e,:) = R\(Q'*b);
                case 4 % Hybrid
                    Htilda = Wsqrt*H;
                    p=colamd(Htilda);
                    R=qr(Htilda(:,p),0);
                    b = Htilda' * Wsqrt * (z-h);
                    dx(p,:) = R'*R\(b(p,:)) ;   
                otherwise
            end;
            time = time + toc;

            if max(abs(dx)) <= eps_tol    
                convergence = 1;
                eps_all = abs(dx);
                it_num = i;
                break;
            else
                % [ Theta; V ] - First index is Bus 1 reference.
                theta(2:size(theta,1)) = theta(2:size(theta,1)) + dx(1:size(theta,1)-1,1);
                V = V + dx(size(theta,1):size(dx,1),1);
            end
        end
        
    case 1
        H_sparse = 1;
        [ H ] = f_measJac_H_v2021( V, theta, Y_bus, topo, ind_meas, N_meas, H_decoupled, H_sparse);
        H_w = H' * W;
        G = H_w * H;

        G_th = G(1:topo.nBus-1,1:topo.nBus-1);
        G_u = G(topo.nBus:2*topo.nBus-1,topo.nBus:2*topo.nBus-1);
        
        L_th = chol(G_th);
        L_u = chol(G_u);
        for i = 1:Max_iter
            [ h ] = f_measFunc_h_v2021( V, theta, Y_bus, topo, ind_meas, N_meas);
            tic;
            rhs = H_w * (z - h);
            rhs_th = rhs(1:topo.nBus-1);
            rhs_u = rhs(topo.nBus:2*topo.nBus-1);
            
            dx_th = L_th\(L_th'\rhs_th);
            dx_u = L_u\(L_u'\rhs_u);
            time = time + toc;
            if max(abs(dx_u)) <= eps_tol && max(abs(dx_th)) <= eps_tol
                convergence = 1;
                it_num = i;
                break;
            else
                % [ Theta; V ] - First index is Bus 1 reference.
                theta(2:size(theta,1)) = theta(2:size(theta,1)) + dx_th;
                V = V + dx_u;
            end
        end        
    otherwise
end
end

