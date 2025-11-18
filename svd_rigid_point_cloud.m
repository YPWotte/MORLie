addpath(genpath('./Auxilliary/ExpMaps'));

N_Array = 50; % [10, 30, 50, 100, 200]; % Number of points
M_Array = 1000; % Number of time-instances per trajectory
L_Array = [10, 30, 50, 100, 200]; % 50; % Number of trajectories

pointcloudplot(50,1000,[10, 30, 50, 100, 200],'Singular values for increasing no. of trajectories');
pointcloudplot([10, 30, 50, 100, 200],1000,50,'Singular values for increasing no. of particles');


function [out] = pointcloudplot(N_Array,M_Array,L_Array,titlestr)
Leg = {};
ind = 1;
figure()
hold on

for N = N_Array
    for M = M_Array
        for L = L_Array
                
            T = 1; % Time
            dt = T/M; % Time step
            
            % Initial point cloud
            mu_pc = [0;0;0];
            var_pc = 1;
            X0 = randn(N*3,L)*var_pc;
            
            % Noise
            var_noise = 0.001; % 0 for no noise at all
            
            
            % Trajectory on SE(3):
            q0 = [rand(3,1);0;0;0]; %rand(6,1);
            q1 = [rand(3,1);0;0;0]; %rand(6,1);
            q = @(t) t*q0 + sin(t)*q1;
            H = @(t) exp_se3(skew(q(t))); % same H for all initial conditions
            
            % Assemble point cloud
            S = zeros(N*3,M*L);
            for k = 1:L
                Sk = zeros(N*3, M);
                Sk(:,1) = X0(:,k); % Initial condition
                for j = 2:M
                    Sk(:,j) = Phi(H(dt*j),Sk(:,1)) + randn(N*3,1)*var_noise;
                end
                S(:,(k-1)*M+1:k*M) = Sk;
            end
            
            [U,SingVals,V] = svd(S,'econ');
            plot(diag(SingVals),'.')
            Leg{ind} = sprintf('Case for %i Particles, %i Trajectories',N,L);
            ind = ind + 1;
        end
    end
end

title(titlestr)
xlabel('Index')
ylabel('Value')
legend(Leg)
hold off
grid on
grid minor
end

%% Plotting
% l = 1; % chosen trajectory
% for j = 1:M 
%     Xj = reshape(S(:,(l-1)*M+j),[3 N]);
%     x = Xj(1,:); y = Xj(2,:); z = Xj(3,:);
%     scatter3(x,y,z);
%     pause(dt);
% end
% 


% Action of SE(3) on N rigidly connected particles
function [Xout] = Phi(H,X)
    [N3,~] = size(X); N = N3/3;
    Xout = H*[reshape(X,[3 N]); ones(1,N)];
    Xout = reshape(Xout(1:3,1:N),[N*3 1]);
end
            