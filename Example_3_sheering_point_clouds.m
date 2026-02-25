%% Algorithm initialized with 2 sheering point clouds, but is able to handle arbitrary amounts.


%% 0. Setup
addpath(genpath('./Auxilliary'));

rng(1);                                                                     % Random seed
NR = 2;                                                                     % Number of rigid bodies
N = 100;                                                                    % Number of points per rigid body
Np = N*NR;                                                                  % Total number of points
M = 1000;                                                                   % Number of time-instances per trajectory
L = 10;                                                                     % Number of trajectories
T = 5;                                                                      % Time
dt = T/M;                                                                   % Time step
mu_pc = [0;0;0];                                                            % Center of first cluster
var_pc = 1;                                                                 % Variation in clusters
var_noise = 1e-5;                                                           % Noise on individual point trajectories 


% 0.1 Initial point clouds
q_offset = [0;0;0;12;0;0];                                                  % Offset between individual point clouds (finite screw)
H_offset = cell(NR,1); H_offset{1} = eye(4);                                 % Offset as matrix in SE(3)
for nr = 2:NR
    H_offset{nr} = H_offset{nr-1}*exp_se3(skew(q_offset));                   
end

X0 = randn(Np*3,L)*var_pc;                                                  % Initial state of point cloud in R^3
for nr = 1:NR                                                               % Applying the offset to separate clusters
    for k = 1:L
        % q_offset = [0;0;0;4;0;0];
        % H_offset = exp(skew(nr*q_offset));
        X0(N*(nr-1)*3+1:N*3*nr,k) = Phi(H_offset{nr}, X0(N*(nr-1)*3+1:N*3*nr,k));
    end
end
% % 0.1.1 Plot initial point-clouds: 
% X = reshape(X0(:,1),[3,Np]);
% x = X(1,:); y = X(2,:); z = X(3,:);
% scatter3(x,y,z);


% 0.2 Trajectories on Aff(3):

% Infinitessimal twist of cluster:
A0 = diag(1:6);
A1 = [1;2;3;4;5;6];
A = @(t) A0*cos(A1*t);                                                      % se(3) component of left twist
B0 = [1,0.5,0.2,0 ; 0.5,2,0.5,0; 0.2,0.5,3,0; 0,0,0,0]/T;                   % sheering and expanding component

% Alternative, completely random twist:
% A0 = rand(6,1); % [rand(3,1); 0;0;0]/T;
% A1 = rand(6,1); % [rand(3,1); 0;0;0]/T; 
% A = @(t) A0 + cos(t)*A1;                                                  
% B0 = [rand(3,4); 0,0,0,0]/T;

% Affine group element H in Aff(3) describing cluster motions
dHdt = @(t,H) reshape(reshape(H,4,4)*(skew(A(t))+B0*sin(2*t)),[],1);        % ODE on GL(4)
Hsol = ode45(dHdt,[0,T],reshape(eye(4),[],1));
H = @(t,nr) H_offset{nr}*reshape(deval(Hsol,t),4,4)*inv_se3(H_offset{nr});  % same H for all initial conditions

% dqdt = @(t,q) K(q,'inverse')*A(t);                                        % Alternatively, ODE in local exponential coordinates (works on SE(3) only)
% qsol = ode45(dqdt,[0,T],zeros(6,1));
% q = @(t) deval(qsol,t);
% H = @(t,nr) H_offset{nr}*exp_se3(skew(q(t)))*inv_se3(H_offset{nr});        

% 0.3 Assemble point cloud trajectories with and without noise
S = zeros(Np*3,M*L);                                                        % Noisy measurement
S_GT = zeros(Np*3,M*L);                                                     % Ground truth, without noise
for k = 1:L                                                                 % Compute noisy measurements and ground truth
    Sk = zeros(Np*3, M);                                                    % Noisy measurement on k-th trajectory
    Sk_GT = zeros(Np*3, M);                                                 % Ground truth on k-th trajectory
    Sk(:,1) = X0(:,k);                                                      % Set initial condition
    Sk_GT(:,1) = X0(:,k);
    for j = 2:M                                                             % Apply affine motion of cluster and add noise
        for nr = 1:NR
            Sk((nr-1)*3*N+1:nr*3*N,j) = Phi(H(dt*(j-1),nr),Sk((nr-1)*3*N+1:nr*3*N,1)) + randn(N*3,1)*var_noise;
            Sk_GT((nr-1)*3*N+1:nr*3*N,j) = Phi(H(dt*(j-1),nr),Sk_GT((nr-1)*3*N+1:nr*3*N,1));
        end
    end
    S(:,(k-1)*M+1:k*M) = Sk; 
    S_GT(:,(k-1)*M+1:k*M) = Sk_GT;
end

% 0.4 Ground truth for twist
rho_GT = zeros(6*NR,(M-1)*L);                                               % infinitessimal right twist at time dt*j
rho_GT_fin = zeros(6*NR,(M-1)*L);                                           % piece-wise constant right finite twist at time dt*j
for k = 1:L
    for j = 1:(M-1)
        for nr = 1:NR
            rho_GT((nr-1)*6+1:nr*6,(k-1)*(M-1)+j) = AdH(H_offset{nr})*AdH(H((j-1)*dt,1))*A((j-1)*dt);
            rho_GT_fin((nr-1)*6+1:nr*6,(k-1)*(M-1)+j) = AdH(H((j-1)*dt,nr))*unskew(log_se3(H((j-1)*dt,nr)\H(j*dt,nr)))/dt;
        end
    end
end


%% 1.0 MORLie: Pre-processing (Algorithm 2)

%% 1.1 Fitting Infinitessimal Generators to Vector fields (Algorithm 2: Velocity-based variant)

S_red = zeros(Np*3,(M-1)*L);                                                % S with final point removed from each trajectory  
Sp = zeros(Np*3,(M-1)*L);                                                   % Difference between successive points (on the same trajectory)
for k = 1:L                                                                 % Loop over all trajectories 
    Sk = S(:,(k-1)*M+1:k*M);                                                
    Sk_red = Sk(:,1:end-1);
    Spk = Sk(:,2:end) - Sk(:,1:end-1);
    S_red(:,(k-1)*(M-1)+1:k*(M-1)) = Sk_red;  
    Sp(:,(k-1)*(M-1)+1:k*(M-1)) = Spk;                                      % Effectively contains dX = X_+ - X;
end

% 1.1.2 Pick random point and its nearest neighbors as local cluster
fprintf("Beginning velocity-based identification ... \n");
Time_FitInfGen = tic;
Nc = 10;                                                                    % Neighbors per cluster
NR_Lim = 10;                                                                % Max number of rigid bodies to identify
NR_Id_Inf = 0;                                                              % Number of identified rigid bodies, updated later
dim_G = 12;                                                                 % Dimension of group
m_G = 4;                                                                    % Dimension of group representation

ixNN_Inf = zeros(L,NR_Lim);                                                 % nearest neighbor callibration index
slack = 0.45;                                                               % Classification slack. Between -0.5 and 0.5. Points that are in a given cluster with probability (0.5-slack) are associated with that cluster.
ctol = sqrt(N)*(var_noise+(var_noise==0))/dt;                                 % Tolerance for assigning points to a cluster

rho_loc = zeros(dim_G,(M-1)*L);                                             % Store generator of local cluster
rho_Data_Inf = zeros(dim_G,(M-1)*L);                                        % Store generators of all clusters
IndOut = true(L,Np);                                                        % Outgroup: points not yet assigned to a cluster
IndCand_Inf = cell(NR_Lim,1);                                               % Candidate clusters
for nr = 1:NR_Lim
    IndCand_Inf{nr} = false(L,Np);                                          
end

for k = 1:L
    notConverged = true;
    nr = 0;
    NIndCand_k = 0;                                                         % Total number of points that are already assigned to candidate clusters, updated later
    IndOut_k = IndOut(k,:);
    while notConverged
        nr = nr+1;
        IndCand_k = IndCand_Inf{nr}(k,:);
        ixNN_Inf(k,nr) = ceil(rand(1)*(Np-NIndCand_k));

        % 1.1.3 Fit infinitessimal generator to velocity of local cluster
        for j = 1:(M-1)
            ind = (k-1)*(M-1)+j;

            % All particles at current time instance:
            X = S_red(:,ind); 
            dX = Sp(:,ind); 

            % Reshape X and dX to identify local cluster
            S_red_kNN = reshape(X,[3, Np]);                                 % All particles
            Sp_red_kNN = reshape(dX,[3, Np]); 
            S_out_kNN = S_red_kNN(:,IndOut_k);                              % Out-group
            Sp_out_kNN = Sp_red_kNN(:,IndOut_k); 
            X_Out = reshape(S_out_kNN, [], 1);                              % Re-shaped out-group
            dX_Out = reshape(Sp_out_kNN, [], 1);

            % Identify local cluster (i.e., nearest neighbors)
            xNN = S_out_kNN(:,ixNN_Inf(k,nr));                              % VIP (very important particle)
            Idx = knnsearch(S_out_kNN.',xNN.',"K",Nc);                      % Indices of neighbor of VIP
            X_loc = reshape(S_out_kNN(:,Idx),[], 1);                        % Neighbors of VIP (including VIP)
            dX_loc = reshape(Sp_out_kNN(:,Idx),[], 1);                      % Velocity of neighbors 

            % Find velocity of local cluster
            A_loc = pinv(JacPhiGL(X_loc))*dX_loc/dt;                        % Local infiniessimal generator
            rho_loc((nr-1)*dim_G+1:nr*dim_G,ind) = A_loc;                   
    
            % 1.1.4 Extend local cluster to candidate rigid body cluster
            diffInfGen = InfGenGL(skewgl(A_loc),X) - dX/dt; 
            error = vecnorm(reshape(diffInfGen,3,[]));                      % Extent to which all particles fail to follow local cluster

            % Identify Candidate rigid body cluster
            IndCand_j = boolean((error < ctol).*(IndOut_k));                % Indices of candidate rigid body cluster, candidates restricted to outgroup
            NCand_j = sum(IndCand_j);                                       % Number of points in candidate rigid body cluster

            % 1.1.5 Compute more accurate infinitessimal generator based on candidate rigid body cluster
            XCand_j = reshape(S_red_kNN(:,IndCand_j),[NCand_j*3,1]);        
            dXCand_j = reshape(Sp_red_kNN(:,IndCand_j),[NCand_j*3,1]);      
            rho_Data_Inf((nr-1)*dim_G+1:nr*dim_G,ind) = ...
                pinv(JacPhiGL(XCand_j))*dXCand_j/dt; 
            
            IndCand_k = IndCand_k + IndCand_j/(M-1);                        % Average association of points with candidate cluster (close to 1 for points that are often in the cluster)
        end
        IndCand_Inf{nr}(k,:) = boolean(round(IndCand_k+slack));             % Binary decision to determine candidate cluster
        IndOut_k = boolean(IndOut_k.*(1-IndCand_Inf{nr}(k,:)));
        NIndCand_k = NIndCand_k + sum(round(IndCand_k+slack));              % Total number of points that are already assigned to candidate clusters
        if nr == NR_Lim                                                     % Failure, too many clusters
            notConverged = false;
            fprintf("Velocity-based fit: failed to converge on trajectory %i, reached cluster limit with %i particles not assigned. \n",k,N - NIndCand_k);
        elseif NIndCand_k == Np                                             % Success, converged
            notConverged = false;
            fprintf("Velocity-based fit: converged on trajectory %i, identified %i clusters. \n",k, nr);
        end
    end
    % 1.1.6 Repeat steps from 1.1.2 for unassigned points, until none remain
    if nr > NR_Id_Inf
        NR_Id_Inf = nr;
    end
end

Time_FitInfGen = toc(Time_FitInfGen);
fprintf("Finished velocity-based fit after T = %i s. \n",Time_FitInfGen)

%% 1.2 Fit finite generator (Algorithm 2: Velocity-free variant, corresponding to article)
fprintf("Beginning velocity-free identification ... \n");
Time_FitFinGen = tic;

rho_loc = zeros(dim_G,(M-1)*L);                                             % Store generator of local cluster
rho_Data_Fin = zeros(dim_G,(M-1)*L);                                        % Store generators of all clusters
ctol = sqrt(N)*(var_noise+(var_noise==0));                                  % Tolerance for assigning points to a cluster
% q_init = zeros(dim_G,1);
% options = optimoptions('fsolve','Display','none','MaxIterations',10,'Algorithm','levenberg-marquardt');

IndCand_Fin = cell(NR_Lim,1);
IndOut = true(L,Np);
for nr = 1:NR_Lim
    IndCand_Fin{nr} = false(L,Np);
end
ixNN_Fin = zeros(L,NR_Lim);                                                 % nearest neighbor callibration index
slack = 0.45;                                                               % Classification slack. Between -0.5 and 0.5. Points that are in a given cluster with probability (0.5-slack) are associated with that cluster.
NR_Id_Fin = 0;                                                              % Number of identified rigid bodies, updated later

nransac = 1;                                                                % number of random initial particles picked per search for local cluster iteration
for k = 1:L
    notConverged = true;
    nr = 0;
    NIndCand_k = 0;                                                         % Total number of points that are already assigned to candidate clusters, updated later
    IndOut_k = IndOut(k,:);
    while notConverged
        nr = nr+1;
        max_identified = 0;
        IndCand_k = false(1,Np);
        IndCand_k_best = IndCand_Fin{nr}(k,:);
        for iransac = 1:nransac
            ixNN_cand = ceil(rand(1)*(Np-NIndCand_k));
            % 1.2.3 Fit finite generator to motion of local cluster
            for j = 1:(M-1)
                ind = (k-1)*(M-1)+j;

                % All particles at current time instance:
                X = S_red(:,ind);                                           
                dX = Sp(:,ind);                                             
                X_next = X + dX;    

                % Reshape X and dX to identify local cluster
                S_red_kNN = reshape(X,[3, Np]);                             % All particles
                Sp_red_kNN = reshape(dX,[3, Np]); 
                S_out_kNN = S_red_kNN(:,IndOut_k);                          % Out-group
                Sp_out_kNN = Sp_red_kNN(:,IndOut_k); 
                X_Out = reshape(S_out_kNN, [], 1);                          % Re-shaped out-group
                dX_Out = reshape(Sp_out_kNN, [], 1);

                % Identify local cluster (i.e., nearest neighbors)
                xNN = S_out_kNN(:,ixNN_cand);
                Idx = knnsearch(S_out_kNN.',xNN.',"K",Nc);
                X_loc = reshape(S_out_kNN(:,Idx),[], 1);
                dX_loc = reshape(Sp_out_kNN(:,Idx),[], 1);

                % Find generator of local cluster motion
                X_next_loc = X_loc+dX_loc;
                A_loc = loggl(LM2_G(@(G,X) Phi(G,X), dim_G, @(A)skewgl(A), ...
                    @(A)expgl(A), eye(4), X_loc, X_next_loc, 2, 0));  
                rho_loc((nr-1)*dim_G+1:nr*dim_G,ind) = A_loc;

                % 1.2.4 Extend local cluster to candidate rigid body cluster
                diffFinGen = Phi(expgl(A_loc),X) - X_next; 
                error = vecnorm(reshape(diffFinGen,3,[]));                  % Extent to which all particles fail to follow local cluster

                % Identify Candidate rigid body cluster
                IndCand_j = boolean((error < ctol).*(IndOut_k));            % Indices of candidate rigid body cluster, candidates restricted to outgroup
                NCand_j = sum(IndCand_j);                                   % Number of points in candidate rigid body cluster

                % 1.2.5 Compute more accurate infinitessimal generator based on candidate rigid body cluster
                XCand_j = reshape(S_red_kNN(:,IndCand_j),[NCand_j*3,1]);
                dXCand_j = reshape(Sp_red_kNN(:,IndCand_j),[NCand_j*3,1]);
                XCand_j_next = XCand_j + dXCand_j;
                rho_Data_Fin((nr-1)*dim_G+1:nr*dim_G,ind) = loggl(LM2_G(@(G,X) Phi(G,X),dim_G,@(A)skewgl(A),@(A)expgl(A),eye(4),XCand_j,XCand_j_next,4,0))/dt;%fsolve(@(q) norm(Phi(expgl(q),XCand_j) - XCand_j_next),q_init,options)/dt;   
                % f = @(G) norm(Phi(G,XCand_j) - XCand_j_next); smin = fminsearch(@(A)f(expgl(A)),zeros(12,1));
                if norm(rho_Data_Fin((nr-1)*dim_G+1:nr*dim_G,ind)) == 0
                    FailBool = true;
                    rho_Data_Fin((nr-1)*dim_G+1:nr*dim_G,ind) = loggl(LM2_G(@(G,X) Phi(G,X),dim_G,@(A)skewgl(A),@(A)expgl(A),eye(4),XCand_j,XCand_j_next,4,0))/dt;%fsolve(@(q) norm(Phi(expgl(q),XCand_j) - XCand_j_next),q_init,options)/dt;   
                end
                

                IndCand_k = IndCand_k + IndCand_j/(M-1);                    % average association of points with candidate cluster (close to 1 for points that are often in the cluster)

            end
            NCand_k = sum(boolean(round(IndCand_k+slack)));
            if NCand_k > max_identified
                max_identified = NCand_k;
                ixNN_Fin(k,nr) = ixNN_cand;
                IndCand_k_best = IndCand_k;
            end
            IndCand_k = false(1,Np);
        end
        IndCand_Fin{nr}(k,:) = boolean(round(IndCand_k_best+slack));        % binary decision to determine candidate cluster
        IndOut_k = boolean(IndOut_k.*(1-IndCand_Fin{nr}(k,:)));
        NIndCand_k = NIndCand_k + sum(boolean(round(IndCand_k_best+slack))); % total number of points that are already assigned to candidate clusters
        if nr == NR_Lim                                                     % Failure, too many clusters
            notConverged = false;
            fprintf("Velocity-free fit: failed to converge on trajectory %i, reached cluster limit with %i particles not assigned. \n",k,N - NIndCand_k);
        elseif NIndCand_k == Np % Success, converged
            notConverged = false;
            fprintf("Velocity-free fit: converged on trajectory %i, identified %i clusters. \n",k, nr);
        end
    end
    % 1.2.6 Repeat steps from 1.2.3 for unassigned points, until none remain
    if nr > NR_Id_Fin
        NR_Id_Fin = nr;
    end
end

Time_FitFinGen = toc(Time_FitFinGen);
fprintf("Finished velocity-free fit after T = %i s. \n",Time_FitFinGen)

%% 1.3 Cleaning up the data: uniquely identify clusters (also part of Algorithm 1)

% Infinitessimal generators
p_av_Inf = cell(NR_Id_Inf,L);                                               % Average position of clusters
Ind_Inf = cell(NR_Id_Inf,1);                                                % Final indices of clusters
rho_Sort_Inf = zeros(size(rho_Data_Inf));                                   % Sorted infinitessimal generators
for k = 1:L                                                                 % Identify average cluster positions
    for nr = 1:NR_Id_Inf
        ind = (k-1)*(M-1)+1;
        X = S_red(:,ind);
        S_red_kNN = reshape(X,[3, Np]);
        IndCand_k_Inf = IndCand_Inf{nr}(k,:);
        p_av_Inf{nr,k} = mean(S_red_kNN(:,IndCand_k_Inf),2);
    end
end

for k = 1:L                                                                 % Group nearest clusters using their average positions
    ind = (k-1)*(M-1)+1;
    indend = k*(M-1);
    for nr_cluster = 1:NR_Id_Inf
        min_dist = inf;
        for nr = 1:NR_Id_Inf
            cluster_dist = norm(p_av_Inf{nr_cluster,k} - p_av_Inf{nr,1});
            if cluster_dist < min_dist
                min_dist = cluster_dist;
                Ind_Inf{nr_cluster}(k,:) = IndCand_Inf{nr}(k,:);
                rho_Sort_Inf((nr_cluster-1)*dim_G+1:nr_cluster*dim_G,ind:indend) = rho_Data_Inf((nr-1)*dim_G+1:nr*dim_G,ind:indend) ;
            end
        end
    end
end


% Finite generators
p_av_Fin = cell(NR_Id_Fin,L);                                               % Average position of clusters
Ind_Fin = cell(NR_Id_Fin,1);                                                % Final indices of clusters
rho_Sort_Fin = zeros(size(rho_Data_Fin));                                   % Sorted infinitessimal generators
for k = 1:L                                                                 % Identify average cluster positions
    for nr = 1:NR_Id_Fin
        ind = (k-1)*(M-1)+1;
        X = S_red(:,ind);
        S_red_kNN = reshape(X,[3, Np]);
        IndCand_k_Fin = IndCand_Fin{nr}(k,:);
        p_av_Fin{nr,k} = mean(S_red_kNN(:,IndCand_k_Fin),2);
    end
end

for k = 1:L                                                                 % Group nearest clusters using their average positions
    ind = (k-1)*(M-1)+1;
    indend = k*(M-1);
    for nr_cluster = 1:NR_Id_Fin
        min_dist = inf;
        for nr = 1:NR_Id_Fin
            cluster_dist = norm(p_av_Fin{nr_cluster,k} - p_av_Fin{nr,1});
            if cluster_dist < min_dist
                min_dist = cluster_dist;
                Ind_Fin{nr_cluster}(k,:) = IndCand_Fin{nr}(k,:);
                rho_Sort_Fin((nr_cluster-1)*dim_G+1:nr_cluster*dim_G,ind:indend) = rho_Data_Fin((nr-1)*dim_G+1:nr*dim_G,ind:indend) ;
            end
        end
    end
end

%% 1.4 SVD: identify subgroups for each cluster (Algorithm 1 in articl)
% 1.4.1: Velocity-based variant
fprintf("Identifying subgroups for velocity-based fit. \n");
Ei_Inf = cell(NR_Id_Inf);
dEi_Inf = cell(NR_Id_Inf);

svd_cutoff_inf = 1e2;                                                       % Percentile cut-off                                                 
decimal_accuracy = 3;
for nr = 1:NR_Id_Inf
    [UInf,SInf,VInf] = svd(rho_Sort_Inf((nr_cluster-1)*dim_G+1:nr_cluster*dim_G,:),"econ");
    diagS = diag(SInf);
    SumS = sum(diagS);
    SumSi = 0; i = 0;
    Indx_subgroup_inf = false(size(diagS));%(max(diagS)/svd_cutoff_inf);
    while (SumSi <= SumS*(1-1/svd_cutoff_inf))
        i = i+1;
        SumSi = SumSi + diagS(i);
        Indx_subgroup_inf(i) = true;
    end
    n = sum(Indx_subgroup_inf);
    Ei_Inf{nr} = cell(n,1);
    for i = 1:n
        Ei_Inf{nr}{i} = skewgl(round(UInf(:,i),decimal_accuracy));
    end
    Ei_Inf{nr} = completeAlgebra(Ei_Inf{nr});
    dEi_Inf{nr} = dual(Ei_Inf{nr});
    fprintf("Cluster %i evolves on %i dimensional subgroup. \n",nr,size(Ei_Inf{nr},1));
%     name_group_nr_inf = sprintf("G%i_inf",nr);
%     makeG(Ei_Inf{nr},name_group_nr_inf);
end

% 1.4.2 Velocity-free variant
fprintf("Identifying subgroups for velocity-free fit. \n");
Ei_Fin = cell(NR_Id_Fin);
dEi_Fin = cell(NR_Id_Fin);
svd_cutoff_fin = 1e2;                                                       % Percentile cut-off
for nr = 1:NR_Id_Inf
    [UFin,SFin,VFin] = svd(rho_Sort_Fin((nr_cluster-1)*dim_G+1:nr_cluster*dim_G,:),"econ");
    % Indx_subgroup_fin = diag(SFin)>(max(diag(SFin))/svd_cutoff_fin);
    diagS = diag(SFin);
    SumS = sum(diagS);
    SumSi = 0; i = 0;
    Indx_subgroup_fin= false(size(diagS));%(max(diagS)/svd_cutoff_inf);
    while (SumSi <= SumS*(1-1/svd_cutoff_fin))
        i = i+1;
        SumSi = SumSi + diagS(i);
        Indx_subgroup_fin(i) = true;
    end
    n = sum(Indx_subgroup_fin);
    Ei_Fin{nr} = cell(n,1);
    for i = 1:n
        Ei_Fin{nr}{i} = skewgl(round(UInf(:,i),decimal_accuracy));
    end
    Ei_Fin{nr} = completeAlgebra(Ei_Fin{nr});
    dEi_Fin{nr} = dual(Ei_Fin{nr});
    fprintf("Cluster %i evolves on %i dimensional subgroup. \n",nr,size(Ei_Fin{nr},1));
% % Optionally generate the group and save it (perfunctory) 
%     name_group_nr_fin = sprintf("G%i_fin",nr);
%     makeG(Ei_Fin{nr},name_group_nr_fin);
end


%% 1.5 Redo fitting for 1D subgroup

% [UInf,SInf,VInf] = svd(rho_Sort_Inf);
% rho_Sort_Inf = UInf(:,1)*(UInf(:,1).'*rho_Sort_Inf);

%% 1.6: Fitting reduced vector field to rho_Data

% 1.4.1 Fit parametric expression to rho_Data
fprintf("Beginning fitting to term-wise minimizers ... \n");
TimeFitting = tic;
T_red = repmat((0:M-2)*dt,1,L);
% % rho_candidate = @(params,t) reshape(params(1:6*NR_Id),[6*NR_Id,1]) + reshape(params((6*NR_Id*3)+1:(6*NR_Id*4)),[6*NR_Id,1])*t + reshape((6*NR_Id*4)+1:(6*NR_Id*5),[6*NR_Id,1])*t.^2 + reshape((6*NR_Id*1)+1:(6*NR_Id*2),[6*NR_Id,1])*sin(t) + reshape((6*NR_Id*2)+1:(6*NR_Id*3),[6*NR_Id,1])*cos(t); % Candidate rho
% rho_candidate_single = @(params,t) reshape(params(1:6),[6,1]) + reshape(params(19:24),[6,1])*t + reshape(params(25:30),[6,1])*t.^2 + reshape(params(7:12),[6,1])*sin(t) + reshape(params(13:18),[6,1])*cos(t); % Candidate rho
% rho_candidate = @(params,t) [rho_candidate_single(params(1:end/2),t); rho_candidate_single(params(end/2+1:end),t)];
% 
% params0 = zeros(30*NR_Id,1);
% 
% [params,resnorm,residual,exitflag,output,lambda,jacobian] =...
%    lsqcurvefit(rho_candidate,params0,T_red,rho_Data);
% TimeFitting = toc(TimeFitting);
% rho_param = @(t) rho_candidate(params,t);

% 1.4.2 Truncated Fourier series to rho_Data
Average = zeros(dim_G*NR_Id_Inf,M-1);
for k = 1:L
    Average = Average + rho_Sort_Inf(:,(k-1)*(M-1)+1:k*(M-1))/L;
end

% 1.4.3 Fit cubic spline to rho_Data
pp = spline((0:M-2)*dt,Average); 
rho_fit_spline = @(t) ppval(pp,t);

% 1.4.3b Fit smoothed cubic spline to data
Average_trunc = [Average(:,1:M/10:end),Average(:,end)];
pp_smoothed = spline((0:M/10:M)*dt,Average_trunc); 
pp_smoothed_coefs = pp_smoothed.coefs;
rho_fit_spline_candidate = @(pp_coefs,t) ppval(coefs_to_pp(pp_smoothed,pp_coefs),t);

[pp_smoothed_coefs,resnorm,residual,exitflag,output,lambda,jacobian] =...
   lsqcurvefit(rho_fit_spline_candidate,pp_smoothed_coefs,T_red,rho_Sort_Inf);
rho_fit_spline_smooth = @(t) rho_fit_spline_candidate(pp_smoothed_coefs,t);

% 1.4.4 Fit Hermite spline to rho_Data
pp_chip = pchip((0:M-2)*dt,Average);
rho_chip = @(t) ppval(pp_chip,t);

% 1.4.4b Fit smoothed Hermite spline to rho_Data
pp_chips = spline((0:M/10:M)*dt,Average_trunc); 
pp_coefs = pp_chips.coefs;
rho_fit_chip_candidate = @(pp_coefs,t) ppval(coefs_to_pp(pp_chips,pp_coefs),t);

[pp_coefs,resnorm,residual,exitflag,output,lambda,jacobian] =...
   lsqcurvefit(rho_fit_chip_candidate,pp_coefs,T_red,rho_Sort_Inf);
rho_chips = @(t) rho_fit_chip_candidate(pp_coefs,t);


% 1.4.5 Truncated Fourier series to rho_Data
Average_Fin = zeros(dim_G*NR_Id_Fin,M-1);
for k = 1:L
    Average_Fin = Average_Fin + rho_Sort_Fin(:,(k-1)*(M-1)+1:k*(M-1))/L;
end

% 1.4.6 Fit cubic spline to rho_Data
ppfin = spline((0:M-2)*dt,Average_Fin); 
rho_fit_spline_fin = @(t) ppval(ppfin,t);

% 1.4.6b Fit smoothed cubic spline to data
Average_trunc_fin = [Average_Fin(:,1:M/10:end),Average_Fin(:,end)];
pp_smoothed_fin = spline((0:M/10:M)*dt,Average_trunc_fin); 
pp_smoothed_coefs_fin = pp_smoothed_fin.coefs;
rho_fit_spline_candidate_fin = @(pp_coefs,t) ppval(coefs_to_pp(pp_smoothed_fin,pp_coefs),t);

[pp_smoothed_coefs_fin,resnorm,residual,exitflag,output,lambda,jacobian] =...
   lsqcurvefit(rho_fit_spline_candidate_fin,pp_smoothed_coefs_fin,T_red,rho_Sort_Fin);
rho_fit_spline_smooth_fin = @(t) rho_fit_spline_candidate_fin(pp_smoothed_coefs_fin,t);
% rho_fit_spline = @(t) ppval(pp,t);

% 1.4.7 Fit Hermite spline to rho_Data
pp_chip_fin = pchip((0:M-2)*dt,Average_Fin);
rho_chip_fin = @(t) ppval(pp_chip_fin,t);

% 1.4.4b Fit smoothed Hermite spline to rho_Data
pp_chips_fin = spline((0:M/10:M)*dt,Average_trunc_fin); 
pp_coefs_fin = pp_chips_fin.coefs;
rho_fit_chip_candidate_fin = @(pp_coefs,t) ppval(coefs_to_pp(pp_chips_fin,pp_coefs),t);

[pp_coefs_fin,resnorm,residual,exitflag,output,lambda,jacobian] =...
   lsqcurvefit(rho_fit_chip_candidate_fin,pp_coefs_fin,T_red,rho_Sort_Fin);
rho_chips_fin = @(t) rho_fit_chip_candidate_fin(pp_coefs_fin,t);

TimeFitting = toc(TimeFitting);
fprintf("Finished fitting term-wise minimizers after T = %i s. \n", TimeFitting);


%% 2.0 Reconstruction: determine Phi: Aff(3) x Aff(3) x R^{3*Np} -> R^{3*Np}

fprintf("Beginning to solve reduced ODEs...\n")
TimeODE = tic;

% 2.1 Various ODEs
dgdt_spline = @(t,g) R_G(g,rho_fit_spline_smooth(t),m_G);  
dgdt_spline_fin = @(t,g) R_G(g,rho_fit_spline_smooth_fin(t),m_G);
dgdt_chip = @(t,g) R_G(g,rho_chips(t),m_G);
dgdt_chip_fin = @(t,g) R_G(g,rho_chips_fin(t),m_G);

% 2.2 Solutions 
g0_Inf = zeros(m_G*m_G*NR_Id_Inf,1);
for nr = 1:NR_Id_Inf
    g0_Inf((nr-1)*m_G*m_G+1:nr*m_G*m_G,:) = reshape(eye(m_G),[],1);
end
sol_spline = ode45(@(t,g) dgdt_spline(t,g), [0 T], g0_Inf);
sol_chip = ode45(@(t,g) dgdt_chip(t,g), [0 T], g0_Inf);

g0_Fin = zeros(m_G*m_G*NR_Id_Fin,1);
for nr = 1:NR_Id_Fin
    g0_Fin((nr-1)*m_G*m_G+1:nr*m_G*m_G,:) = reshape(eye(m_G),[],1);
end
sol_spline_fin = ode45(@(t,g) dgdt_spline_fin(t,g), [0 T], g0_Fin);
sol_chip_fin = ode45(@(t,g) dgdt_chip_fin(t,g), [0 T], g0_Fin);

TimeODE = toc(TimeODE);
fprintf("Solved reduced ODEs after T = %i s. \n", TimeODE);

G_spline = @(t) gvec_to_g(deval(sol_spline,t),m_G);
G_spline_fin = @(t) gvec_to_g(deval(sol_spline_fin,t),m_G);
G_chip = @(t) gvec_to_g(deval(sol_chip,t),m_G);
G_chip_fin = @(t) gvec_to_g(deval(sol_chip_fin,t),m_G);

% 2.3 Reconstructed trajectories
fprintf("Beginning to reconstruct full solutions and compute errors...\n")
TimeReconstruct = tic;

l = 1;                                                                      % chosen trajectory
k = l;
Ind_k_inf = Ind_Inf{k}(k,:);
for nr = 1:NR_Id_Inf
    Ind_k_inf(nr,:) = Ind_Inf{nr}(k,:);
end

Sk_spline = zeros(Np*3, M);
Sk_spline(:,1) = X0(:,k); % Initial condition
for j = 2:M
    Sk_spline(:,j) = PhiG(G_spline(dt*(j-1)),Sk_spline(:,1),Ind_k_inf);
end

Sk_chip = zeros(Np*3, M);
Sk_chip(:,1) = X0(:,k); % Initial condition
for j = 2:M
    Sk_chip(:,j) = PhiG(G_chip(dt*(j-1)),Sk_chip(:,1),Ind_k_inf);
end


Ind_k_fin = Ind_Fin{k}(k,:);
for nr = 1:NR_Id_Inf
    Ind_k_fin(nr,:) = Ind_Fin{nr}(k,:);
end

Sk_spline_fin = zeros(Np*3, M);
Sk_spline_fin (:,1) = X0(:,k); % Initial condition
for j = 2:M
    Sk_spline_fin (:,j) = PhiG(G_spline_fin(dt*(j-1)),Sk_spline_fin(:,1),Ind_k_fin);
end

Sk_chip_fin = zeros(Np*3, M);
Sk_chip_fin(:,1) = X0(:,k); % Initial condition
for j = 2:M
    Sk_chip_fin (:,j) = PhiG(G_chip_fin(dt*(j-1)),Sk_chip_fin(:,1),Ind_k_fin);
end

% 2.4 Errors: average particle distance 
p_norm = 2;
% 2.4.1 Over whole trajectory: 
Sk_GT = S_GT(:,(k-1)*M+1:k*M);
% Ett_param = vecnorm(Sk_param - Sk_GT,norm);
Ett_spline = dist(Sk_spline,Sk_GT,p_norm); 
Ett_chip = dist(Sk_chip,Sk_GT,p_norm); 
Ett_spline_fin = dist(Sk_spline_fin,Sk_GT,p_norm);
Ett_chip_fin = dist(Sk_chip_fin,Sk_GT,p_norm); 

% 2.4.2 Between successive points:
% Esp_param = zeros(M-2,1);
Esp_spline = zeros(M-2,1);
Esp_chip = zeros(M-2,1);
Esp_spline_fin = zeros(M-2,1);
Esp_chip_fin = zeros(M-2,1);
for j = 2:M
    % Esp_param(j-1) = vecnorm(PhiG(expgl(dt*rho_param((j-2)*dt)),Sk_GT(:,j-1),IndCand_k) - Sk_GT(:,j),norm);
    Esp_spline(j-1) = dist(PhiG(expG(dt*rho_fit_spline_smooth((j-2)*dt),dim_G),Sk_GT(:,j-1),Ind_k_inf), Sk_GT(:,j),p_norm); 
    Esp_chip(j-1) = dist(PhiG(expG(dt*rho_chips((j-2)*dt),dim_G),Sk_GT(:,j-1),Ind_k_inf), Sk_GT(:,j),p_norm); 
    Esp_spline_fin(j-1) = dist(PhiG(expG(dt*rho_fit_spline_smooth_fin((j-2)*dt),dim_G),Sk_GT(:,j-1),Ind_k_fin), Sk_GT(:,j),p_norm); 
    Esp_chip_fin(j-1) = dist(PhiG(expG(dt*rho_chips_fin((j-2)*dt),dim_G),Sk_GT(:,j-1),Ind_k_fin), Sk_GT(:,j),p_norm); 
end

TimeReconstruct = toc(TimeReconstruct);
fprintf("Reconstructed for trajectory %i and computed errors after T = %i s.\n",l, TimeReconstruct);


%% 3.0 Plotting
day = datetime('now','Format','y-M-d');
time = datetime('now','Format','H-m-s');
FolderName = sprintf('PseudoRigidClouds_%s_at_%s',day,time);
mkdir(FolderName);

fig = figure();
clear Frames;
skip = 10;
Frames(M/skip) = struct('cdata', [], 'colormap', []);

video_filename = sprintf('ReconstructionVideo_%s_at_%s',day,time);
video_filename = strcat(FolderName,'/',video_filename);
Sl = reshape(S(:,(l-1)*M+1:l*M),[3, Np*M]);
xm = min(Sl(1,:)); xp = max(Sl(1,:));
ym = min(Sl(2,:)); yp = max(Sl(2,:));
zm = min(Sl(3,:)); zp = max(Sl(3,:));
pad = abs(xp-xm)*0.05;
% title('Reconstruction of rigid point cloud')
% legend('Full $\gamma(t)$','Reconstructed $\bar{\gamma} = \Phi\big(g(t),\gamma_0\big)$',interpreter='latex')
% legend('')
ind = 0;
for j = 1:skip:M
    ind = ind+1;
    Xj_GT = reshape(S_GT(:,(l-1)*M+j),[3 Np]);
    Xj_fit = reshape(Sk_chip_fin(:,j),[3 Np]);
    x_GT = Xj_GT(1,:); y_GT = Xj_GT(2,:); z_GT = Xj_GT(3,:);
    x_fit = Xj_fit(1,:); y_fit = Xj_fit(2,:); z_fit = Xj_fit(3,:);
    scatter3(x_GT,y_GT,z_GT,'o');
    hold on
    scatter3(x_fit,y_fit,z_fit,'*');
    hold off
    title('Reconstruction of rigid point cloud')
    legend('Full $P_i(t)$','Reconstructed $\bar{P}_i(t) = \Phi\big(g(t), P_{i,0}\big)$',interpreter='latex',location='northeast')
    xlim([xm-pad xp+pad]); ylim([ym-pad yp+pad]); zlim([zm-pad zp+pad]);
    Frames(ind) = getframe(fig);
end
video = VideoWriter(sprintf('%s.avi',video_filename));
video.FrameRate = 1/dt/skip;
open(video)
writeVideo(video,Frames(2:end));
close(video)


% Plot a few frames:

ind = 0;
for j = [1, M/4, M/2, 3*M/4, M]
    ind = ind+1;
    fig = figure();
    name = sprintf('ReconstructionSnapshot_%i_on_%s_at_%s.svg',ind,day,time);
    name = strcat(FolderName,'/',name);
    Xj_GT = reshape(S_GT(:,(l-1)*M+j),[3 Np]);
    Xj_fit = reshape(Sk_chip_fin(:,j),[3 Np]);
    x_GT = Xj_GT(1,:); y_GT = Xj_GT(2,:); z_GT = Xj_GT(3,:);
    x_fit = Xj_fit(1,:); y_fit = Xj_fit(2,:); z_fit = Xj_fit(3,:);
    scatter3(x_GT,y_GT,z_GT,'o');
    hold on
    scatter3(x_fit,y_fit,z_fit,'*');
    hold off
    title(sprintf('T = %.2f s', j*dt),interpreter="latex");
    lgd = legend('Full $P_i(t)$','Reconstructed $\bar{P}_i(t) = \Phi\big(g(t), P_{i,0}\big)$',interpreter='latex',location='northeast');
    fontsize(17,"points"); % fontsize(lgd,14,'points');
    xlim([xm-pad xp+pad]); ylim([ym-pad yp+pad]); zlim([zm-pad zp+pad]);
    saveas(fig,name);
end

% fig = figure();
% clear Frames;
% Frames(M) = struct('cdata', [], 'colormap', []);
% day = datetime('now','Format','y-M-d');
% time = datetime('now','Format','H-m-s');
% video_filename = sprintf('Rigid_point_cloud_%s_at_%s',day,time);
% Sl = reshape(S(:,(l-1)*M+1:l*M),[3, N*M]);
% xm = min(Sl(1,:)); xp = max(Sl(1,:));
% ym = min(Sl(2,:)); yp = max(Sl(2,:));
% zm = min(Sl(3,:)); zp = max(Sl(3,:));
% pad = abs(xp-xm)*0.05;
% % title('Reconstruction of rigid point cloud')
% % legend('Full $\gamma(t)$','Reconstructed $\bar{\gamma} = \Phi\big(g(t),\gamma_0\big)$',interpreter='latex')
% % legend('')
% for j = 1:M 
%     Xj_GT = reshape(S_GT(:,(l-1)*M+j),[3 N]);
%     Xj_fit = reshape(Sk_chip(:,j),[3 N]);
%     x_GT = Xj_GT(1,:); y_GT = Xj_GT(2,:); z_GT = Xj_GT(3,:);
%     x_fit = Xj_fit(1,:); y_fit = Xj_fit(2,:); z_fit = Xj_fit(3,:);
%     scatter3(x_GT,y_GT,z_GT,'o');
%     hold on
%     scatter3(x_fit,y_fit,z_fit,'*');
%     hold off
%     title('Reconstruction of rigid point cloud')
%     legend('Full $\gamma(t)$','Reconstructed $\bar{\gamma} = \Phi\big(g(t),\gamma_0\big)$',interpreter='latex',location='northeast')
%     xlim([xm-pad xp+pad]); ylim([ym-pad yp+pad]); zlim([zm-pad zp+pad]);
%     Frames(j) = getframe(fig);
% end
% video = VideoWriter(sprintf('%s.avi',video_filename));
% video.FrameRate = 1/dt;
% open(video)
% writeVideo(video,Frames(2:end));
% close(video)


name = sprintf('SVD_%s_at_%s.fig',day,time);
name = strcat(FolderName,'/',name);
name2 = sprintf('SVD_%s_at_%s.svg',day,time);
name2 = strcat(FolderName,'/',name2);
fig = figure();
[U,SingVals,V] = svd(S,'econ');
semilogy(diag(SingVals),'.')
title('Singular values of $S$')
xlabel('Index')
ylabel('Value')
lgd = legend('Singular values'); %legend(sprintf('Case for %i Particles, %i Trajectories',N,L));
fontsize(17,"points"); % fontsize(lgd,14,'points')
grid on 
grid minor
saveas(fig,name);
saveas(fig,name2);

name = sprintf('SVD_VB_%s_at_%s.fig',day,time);
name = strcat(FolderName,'/',name);
name2 = sprintf('SVD_VB_%s_at_%s.svg',day,time);
name2 = strcat(FolderName,'/',name2);
fig = figure();
[U1,SingVals1,V1] = svd(rho_Sort_Inf,'econ');
[U2,SingVals2,V2] = svd(rho_Sort_Fin,'econ');
semilogy(diag(SingVals1),'o')
hold on
semilogy(diag(SingVals2),'o')
hold off
title('Singular values of $S_{aff(3)}$')
xlabel('Index')
ylabel('Value')
lgd = legend('Velocity-based','Velocity-free'); %legend(sprintf('Case for %i Particles, %i Trajectories',N,L));
fontsize(17,"points"); % fontsize(lgd,14,'points'); 
grid on 
grid minor
saveas(fig,name);
saveas(fig,name2);

% name = sprintf('SVD_VF_%s_at_%s.fig',day,time);
% name = strcat(FolderName,'/',name);
% name2 = sprintf('SVD_VF_%s_at_%s.svg',day,time);
% name2 = strcat(FolderName,'/',name2);
% fig = figure();
% [U1,SingVals1,V1] = svd(rho_Sort_Fin,'econ');
% semilogy(diag(SingVals1),'o')
% title('Singular values of velocity-free reduced snapshots')
% xlabel('Index')
% ylabel('Value')
% lgd = legend(sprintf('Case for %i Particles, %i Trajectories',N,L));
% fontsize(17,"points"); % fontsize(lgd,14,'points');
% grid on 
% grid minor
% saveas(fig,name);
% saveas(fig,name2);

name = sprintf('TT_Error_%s_at_%s.svg',day,time);
name = strcat(FolderName,'/',name);
fig = figure();
T_tt = (0:M-1)*dt;
%plot(T_tt,Ett_param,T_tt,Ett_fft,T_tt,Ett_spline,T_tt,Ett_chip,T_tt,Ett_param_fin,T_tt,Ett_fft_fin,T_tt,Ett_spline_fin,T_tt,Ett_chip_fin,'LineWidth',2)
%plot(T_tt,Ett_spline,T_tt,Ett_chip,T_tt,Ett_spline_fin,T_tt,Ett_chip_fin,'LineWidth',2)
plot(T_tt,Ett_chip,T_tt,Ett_chip_fin,'LineWidth',2)


title('Full trajectory error')
xlabel('Time in s')
ylabel('Error')
%legend('Parametric Fit','Fourier Fit', 'Spline Fit','Hermite Fit', 'Parametric Fit: Finite twist', 'Fourier Fit: Finite twist', 'Spline Fit: Finite twist','Hermite Fit: Finite twist')
%legend('Spline Fit','Hermite Fit','Spline Fit Finite','Hermite Fit Finite','Location','northwest')
lgd = legend('Velocity-based','Velocity-free','Location','northwest');
fontsize(17,"points"); % fontsize(lgd,14,'points');

grid on 
grid minor
saveas(fig,name);

name = sprintf('SP_Error_%s_at_%s.svg',day,time);
name = strcat(FolderName,'/',name);
fig = figure();
T_sp = (0:M-2)*dt;
% plot(T_sp,Esp_param,T_sp,Esp_fft,T_sp,Esp_spline,T_sp,Esp_chip,T_sp,Esp_param_fin,T_sp,Esp_fft_fin,T_sp,Esp_spline_fin,T_sp,Esp_chip_fin,'LineWidth',2)
% plot(T_sp,Esp_spline,T_sp,Esp_chip,T_sp,Esp_spline_fin,T_sp,Esp_chip_fin,'LineWidth',2)
plot(T_sp,Esp_chip,T_sp,Esp_chip_fin,'LineWidth',2)

title('Step ahead error')
xlabel('Time in s')
ylabel('Error')
% legend('Parametric Fit','Fourier Fit', 'Spline Fit','Hermite Fit', 'Parametric Fit: Finite twist', 'Fourier Fit: Finite twist', 'Spline Fit: Finite twist','Hermite Fit: Finite twist')
% legend('Spline Fit','Hermite Fit','Spline Fit Finite','Hermite Fit Finite','Location','northwest')
lgd = legend('Velocity-based','Velocity-free','Location','northwest');
fontsize(17,"points"); % fontsize(lgd,14,'points');

grid on 
grid minor
saveas(fig,name);
% 
% % % Noisy video
% fig = figure();
% clear Frames;
% skip = 10;
% Frames(M/skip) = struct('cdata', [], 'colormap', []);
% day = datetime('now','Format','y-M-d');
% time = datetime('now','Format','H-m-s');
% video_filename = sprintf('Noise_2_D_nonintegrable_Pseudo_rigid_point_clouds_%s_at_%s',day,time);
% video_filename = strcat(FolderName,'/',video_filename);
% Sl = reshape(S(:,(l-1)*M+1:l*M),[3, Np*M]);
% xm = min(Sl(1,:)); xp = max(Sl(1,:));
% ym = min(Sl(2,:)); yp = max(Sl(2,:));
% zm = min(Sl(3,:)); zp = max(Sl(3,:));
% pad = abs(xp-xm)*0.05;
% % title('Reconstruction of rigid point cloud')
% % legend('Full $\gamma(t)$','Reconstructed $\bar{\gamma} = \Phi\big(g(t),\gamma_0\big)$',interpreter='latex')
% % legend('')
% ind = 0;
% for j = 1:skip:M
%     ind = ind+1;
%     Xj_GT = reshape(S(:,(l-1)*M+j),[3 Np]);
%     Xj_fit = reshape(Sk_chip(:,j),[3 Np]);
%     x_GT = Xj_GT(1,:); y_GT = Xj_GT(2,:); z_GT = Xj_GT(3,:);
%     x_fit = Xj_fit(1,:); y_fit = Xj_fit(2,:); z_fit = Xj_fit(3,:);
%     scatter3(x_GT,y_GT,z_GT,'o');
%     title('Noisy point cloud measurement')
%     legend('Full $\gamma(t)$',interpreter='latex',location='northeast')
%     xlim([xm-pad xp+pad]); ylim([ym-pad yp+pad]); zlim([zm-pad zp+pad]);
%     Frames(ind) = getframe(fig);
% end
% video = VideoWriter(sprintf('%s.avi',video_filename));
% video.FrameRate = 1/dt/skip;
% open(video)
% writeVideo(video,Frames(2:end));
% close(video)




%% Functions

%% Phi: SE3 x R^3N -> R^3N

% Action of SE(3) on N rigidly connected particles
function [Xout] = Phi(H,X)
    [N3,~] = size(X); N = N3/3;
    Xout = H*[reshape(X,[3 N]); ones(1,N)];
    Xout = reshape(Xout(1:3,1:N),[N*3 1]);
end

% Infinitessimal generator of se(3) on N rigidly connected particles
function [X_A] = InfGen(A,X)
    [N3,~] = size(X); N = N3/3;
    X_A = A*[reshape(X,[3 N]); ones(1,N)];
    X_A = reshape(X_A(1:3,1:N),[N*3 1]);
end

% Jacobian mapping components of A to X_A (dimension 3*N by 6)
function [PhiJac] = JacPhi(X)
    [N3,~] = size(X); N = N3/3;
    Ai = eye(6);
    PhiJac = zeros(N*3,6);
    for i = 1:6
        X_A = skew(Ai(:,i))*[reshape(X,[3 N]); ones(1,N)];
        PhiJac(:,i) = reshape(X_A(1:3,1:N),[N*3 1]);
    end
end

function [error] = dist_Phi(H,X,Xp)
    error = norm(Xp - Phi(H,X));
end

function [error] = dist_Phi_vec(H,X,Xp)
    error = (Xp - Phi(H,X)).^2;
end

%% Phi: G x R^3N -> R^3N
% for G = Aff3 x ... x Aff3

function [G] = expG(A,dim_G)
    [N,~] = size(A); NG = N/dim_G;
    G = cell(NG,1);
    for ng = 1:NG
        ind = (ng-1)*dim_G+1:ng*dim_G;
        G{ng} = expgl(A(ind));
    end
end

function [Out] = invG(G)
    [NG,~] = size(G);
    Out = cell(NG,1);
    for ng = 1:NG
        Out{ng} = inv(G{ng});
    end
end

% function [T] = Kg(q,dqdt,type)
%     [N6,~] = size(q);
%     T = zeros(size(q));
%     for ng = 1:N6/6
%         ind = (ng-1)*6+1:ng*6;
%         T(ind) = K(q(ind),type)*dqdt(ind);
%     end
% end

function [Tout] = AdG(G,T)
    [N6,~] = size(T);
    Tout = zeros(size(T));
    NG = N6/6;
    for ng = 1:NG
        ind = (ng-1)*6+1:ng*6;
        Tout(ind) = unskewgl(G{ng}*skewgl(T(ind))/G{ng});
    end
end

function [Xout] = PhiG(g,X,Ind) % Assumes g a (NG,1) cell of SE(3) matrices, Ind a NG x dim(X) boolean matrix that indicates on which elements of X g{ng} will act
    [N3,~] = size(X); N = N3/3;
    Xr = [reshape(X,[3 N]); ones(1,N)];
    Xout = [zeros(3,N); ones(1,N)];

    [NG,~] = size(g); 
    for ng = 1:NG
        Xout(:,Ind(ng,:)) = g{ng}*Xr(:,Ind(ng,:));
    end
    Xout = reshape(Xout(1:3,1:N),[N*3 1]);
end

function [dgdt] = R_G(gvec,A,m_G)
    Ng = length(gvec);
    NA = length(A);
    nId = Ng/m_G^2; N = NA/nId;
    dgdt = zeros(Ng,1);
    for i = 1:nId
        indA = (i-1)*N+1;
        indendA = i*N;
        indg = (i-1)*Ng/nId+1;
        indendg = i*Ng/nId;
        gi = reshape(gvec(indg:indendg),m_G,m_G);
        Ai = skewgl(A(indA:indendA));
        dgdt(indg:indendg) = reshape(Ai*gi,[],1);
    end
end

function [dgdt] = L_G(gvec,A,m_G)
    Ng = length(gvec);
    NA = length(A);
    nId = Ng/m_G^2; N = NA/nId;
    dgdt = zeros(Ng,1);
    for i = 1:nId
        indA = (i-1)*N+1;
        indendA = i*N;
        indg = (i-1)*Ng/nId+1;
        indendg = i*Ng/nId;
        gi = reshape(gvec(indg:indendg),m_G,m_G);
        Ai = skewgl(A(indA:indendA));
        dgdt(indg:indendg) = reshape(gi*Ai,[],1);
    end
end

function g = gvec_to_g(gvec,m_G)
    N = m_G^2;
    nId = length(gvec)/N;
    g = cell(nId,1);
    for i = 1:nId
        ind = (i-1)*N+1;
        indend = i*N;
        g{i} = reshape(gvec(ind:indend),m_G,m_G);
    end
end

% Infinitessimal generator of aff(3) on N rigidly connected particles
function [X_A] = InfGenG(A,X,Ind,nId)
    [N3,~] = size(X); N = N3/3;
    Xr = [reshape(X,[3 N]); ones(1,N)];
    X_A = zeros(4,N);

    [N,~] = size(A); NG = N/nId;
    for ng = 1:NG
        ind = (ng-1)*12+1:ng*12;
        X_A(:,Ind(ng,:)) = skewgl(A(ind))*Xr(:,Ind(ng,:));
    end
    X_A = reshape(X_A(1:3,1:N),[N*3 1]);
end

% Jacobian mapping components of A to X_A (dimension 3*N by 6)
function [PhiJac] = JacPhiG(X,Ind)
    [N3,~] = size(X); N = N3/3;
    [NG,~] = size(Ind);
    
    Ai = eye(NG*12);
    PhiJac = zeros(N*3,6*NG);
    for i = 1:6*NG
        X_Ai = InfGenG(Ai(:,ng),X,Ind);
        PhiJac(:,i) = X_Ai;
    end
end

function [error] = dist_PhiG(g,X,Xp)
    error = norm(Xp - PhiG(g,X));
end

function [error] = dist_PhiG_vec(g,X,Xp)
    error = (Xp - PhiG(g,X)).^2;
end

% Infinitessimal generator of gl(3) x R³ on N rigidly connected particles
function [X_A] = InfGenGL(A,X)
    [N3,~] = size(X); N = N3/3;
    X_A = A*[reshape(X,[3 N]); ones(1,N)];
    X_A = reshape(X_A(1:3,1:N),[N*3 1]);
end

% Jacobian mapping components of A to X_A (dimension 3*N by 12)
function [PhiJac] = JacPhiGL(X)
    [N3,~] = size(X); N = N3/3;
    Ai = eye(12);
    PhiJac = zeros(N*3,12);
    for i = 1:12
        X_A = skewgl(Ai(:,i))*[reshape(X,[3 N]); ones(1,N)];
        PhiJac(:,i) = reshape(X_A(1:3,1:N),[N*3 1]);
    end
end

function [Atilde] = skewgl(A)
    Atilde = [reshape(A,[3, 4]); [0, 0, 0, 0]];
end

function [A] = unskewgl(Atilde)
    A = reshape(Atilde(1:3,:),[12,1]);
end

function [G] = expgl(A) 
    %G = eye(4) + skewgl(A); % Bad exp, only locally valid (small A)
    G = expm(skewgl(A));
end

function [A] = loggl(G) 
    %A = unskewgl(G - eye(4)); % Bad log, only locally valid (small A)
    A = unskewgl(logm(G));
end

function [Aout] = AdGL(g,A)
    Aout = unskewgl(g*skewgl(A)/g);
end

%% Misc

function [pp] = coefs_to_pp(pp,coefs)
    pp.coefs = coefs;
end