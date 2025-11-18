
%% 0. Setup
addpath(genpath('./Auxilliary'));

rng(1);                                                                     % Random seed
N = 100;                                                                    % Number of points
M = 1000;                                                                   % Number of time-instances per trajectory
L = 10;                                                                     % Number of trajectories
T = 1;                                                                      % Time
dt = T/M;                                                                   % Time step
mu_pc = [0;0;0];                                                            % Center of first cluster
var_pc = 1;                                                                 % Variation in clusters
X0 = randn(N*3,L)*var_pc;                                                   % Initial point cloud
var_noise = 1e-5;                                                           % Noise on individual point trajectories 

% 0.1 Trajectory on SE(3):
A0 = diag(1:6);
A1 = [1;2;3;4;5;6];
A = @(t) A0*cos(A1*t);                                                     
B0 = zeros(4,4);
dqdt = @(t,q) K(q,'inverse')*A(t);

% Alternative, specifying trajectory in finite screws:
% q0 = rand(6,1); %[rand(3,1);0;0;0];
% q1 = rand(6,1); %[rand(3,1);0;0;0];
% q = @(t) t*q0 + sin(t)*q1;
% dqdt = @(t) q0 + cos(t)*q1;
% A = @(t) K(q(t))*dqdt(t); 

qsol = ode45(dqdt,[0,T],zeros(6,1));
q = @(t) deval(qsol,t);
H = @(t) exp_se3(skew(q(t))); 

% 0.4 Assemble point cloud
S = zeros(N*3,M*L);                                                         % Noisy measurement
S_GT = zeros(N*3,M*L);                                                      % Ground Truth
for k = 1:L
    Sk = zeros(N*3, M);
    Sk_GT = zeros(N*3, M);
    Sk(:,1) = X0(:,k);                                                      % Initial condition
    Sk_GT(:,1) = X0(:,k);
    for j = 2:M
        Sk(:,j) = Phi(H(dt*(j-1)),Sk(:,1)) + randn(N*3,1)*var_noise;
        Sk_GT(:,j) = Phi(H(dt*(j-1)),Sk_GT(:,1));
    end
    S(:,(k-1)*M+1:k*M) = Sk;
    S_GT(:,(k-1)*M+1:k*M) = Sk_GT;
end

% 0.5 Ground truth 
% A = @(t) K(q(t))*dqdt(t);
rho_GT = zeros(6,(M-1)*L);                                                  % twist at time dt*j
rho_GT_fin = zeros(6,(M-1)*L);                                              % piece-wise constant finite twist at time dt*j
for k = 1:L
    for j = 1:(M-1)
        rho_GT(:,(k-1)*(M-1)+j) = AdH(H((j-1)*dt))*A((j-1)*dt);
        rho_GT_fin(:,(k-1)*(M-1)+j) = AdH(H((j-1)*dt))*unskew(log_se3(H((j-1)*dt)\H(j*dt)))/dt;
    end
end


%% 1.0 MorLie
% 1.1 Pre-processing:
% 1.1.1 Fitting Infinitessimal Generators to Vector fields

S_red = zeros(N*3,(M-1)*L);                                                 % remove final point 
Sp = zeros(N*3,(M-1)*L);                                                    % difference between successive points (on the same trajectory)
for k = 1:L
    Sk = S(:,(k-1)*M+1:k*M);
    Sk_red = Sk(:,1:end-1);
    Spk = Sk(:,2:end) - Sk(:,1:end-1);
    S_red(:,(k-1)*(M-1)+1:k*(M-1)) = Sk_red;
    Sp(:,(k-1)*(M-1)+1:k*(M-1)) = Spk;
end

Time_FitInfGen = tic;
rho_Data = zeros(6,(M-1)*L);
for k = 1:L
    for j = 1:(M-1)
        ind = (k-1)*(M-1)+j;
        X = S_red(:,ind);
        dX = Sp(:,ind);
        rho_Data(:,ind) = pinv(JacPhi(X))*dX/dt; %min(dist(Phi(H,X),X))
    end
end
Time_FitInfGen = toc(Time_FitInfGen);

% 1.1.2 Fitting Finite Twists to Point Cloud Differences

Time_FitFinGen = tic;
rho_Data_Fin = zeros(6,(M-1)*L);
q_init = zeros(6,1);
options = optimoptions('fsolve','Display','none','MaxIterations',10,'Algorithm','levenberg-marquardt');
for k = 1:L
    for j = 1:(M-1)
        ind = (k-1)*(M-1)+j;
        X = S_red(:,ind);
        dX = Sp(:,ind);
        X_next = X+dX; 
        rho_Data_Fin(:,ind) = unskew(log_se3(LM2_SE3(@(H,X) Phi(H,X),eye(4),X,X_next,2,0))); %fsolve(@(q)dist_Phi(exp_se3(skew(q)),X,X_next),q_init,options); %fminsearch(@(q)dist_Phi(exp_se3(skew(q)), X, X_next),q0); 
        q_init = rho_Data_Fin(:,ind);
    end
end
rho_Data_Fin = rho_Data_Fin/dt;
Time_FitFinGen = toc(Time_FitFinGen);

%% 1.2 Fitting reduced vector field to rho_Data

% 1.2.1 Fit parametric expression to rho_Data
TimeFitting = tic;
T_red = repmat((0:M-2)*dt,1,L);
rho_candidate = @(params,t) reshape(params(1:6),[6,1]) + reshape(params(19:24),[6,1])*t + reshape(params(25:30),[6,1])*t.^2 + reshape(params(7:12),[6,1])*sin(t) + reshape(params(13:18),[6,1])*cos(t); % Candidate rho
params0 = zeros(30,1);

[params,resnorm,residual,exitflag,output,lambda,jacobian] =...
   lsqcurvefit(rho_candidate,params0,T_red,rho_Data);
TimeFitting = toc(TimeFitting);
rho_param = @(t) rho_candidate(params,t);

% 1.2.2 Truncated Fourier series to rho_Data
Average = zeros(6,M-1);
for k = 1:L
    Average = Average + rho_Data(:,(k-1)*(M-1)+1:k*(M-1))/L;
end
SumFFT = fft(Average.').';
[sx,sy] = size(Average);
if mod(sy,2) == 0
    SumFFT_Truncated = [SumFFT(:,1:10),SumFFT(:,end-10:end)];
else
    SumFFT_Truncated = [SumFFT(:,1:10),SumFFT(:,end-8:end)];
end
rho_fit_FFT = @(t) real(FFT_Interpolate1D(t,T,SumFFT));

% 1.2.3 Fit cubic spline to rho_Data
pp = spline((0:M-2)*dt,Average); 
rho_fit_spline = @(t) ppval(pp,t);

% 1.2.3b Fit smoothed cubic spline to data
Average_trunc = [Average(:,1:M/10:end),Average(:,end)];
pp_smoothed = spline((0:M/10:M)*dt,Average_trunc); 
pp_smoothed_coefs = pp_smoothed.coefs;
rho_fit_spline_candidate = @(pp_coefs,t) ppval(coefs_to_pp(pp_smoothed,pp_coefs),t);

[pp_smoothed_coefs,resnorm,residual,exitflag,output,lambda,jacobian] =...
   lsqcurvefit(rho_fit_spline_candidate,pp_smoothed_coefs,T_red,rho_Data);
rho_fit_spline_smooth = @(t) rho_fit_spline_candidate(pp_smoothed_coefs,t);
% rho_fit_spline = @(t) ppval(pp,t);

% 1.2.4 Fit Hermite spline to rho_Data
pp_chip = pchip((0:M-2)*dt,Average);
rho_chip = @(t) ppval(pp_chip,t);

% 1.2.5 Fit parametric expression to rho_Data_Fin
rho_candidate_Fin = @(params,t) reshape(params(1:6),[6,1]) + reshape(params(19:24),[6,1])*t + reshape(params(25:30),[6,1])*t.^2 + reshape(params(7:12),[6,1])*sin(t) + reshape(params(13:18),[6,1])*cos(t); % Candidate rho
params0 = zeros(30,1);

[params_Fin,resnorm,residual,exitflag,output,lambda,jacobian] =...
   lsqcurvefit(rho_candidate,params0,T_red,rho_Data_Fin);
rho_param_Fin = @(t) rho_candidate(params_Fin,t);

% 1.2.6 Truncated Fourier series to rho_Data_Fin
Average_FIN = zeros(6,M-1);
for k = 1:L
    Average_FIN = Average_FIN + rho_Data_Fin(:,(k-1)*(M-1)+1:k*(M-1))/L;
end
SumFFT_FIN = fft(Average_FIN.').';
[sx,sy] = size(Average_FIN);
if mod(sy,2) == 0
    SumFFT_Truncated_FIN = [SumFFT_FIN(:,1:10),SumFFT_FIN(:,end-10:end)];
else
    SumFFT_Truncated_FIN = [SumFFT_FIN(:,1:10),SumFFT_FIN(:,end-8:end)];
end
rho_fit_FFT_Fin = @(t) real(FFT_Interpolate1D(t,T,SumFFT_FIN));

% 1.2.7 Fit cubic spline to rho_Data_Fin
pp_fin = spline((0:M-2)*dt,Average_FIN);
rho_fit_spline_fin = @(t) ppval(pp_fin,t);

% 1.2.8 Fit Hermite spline to rho_Data_Fin
pp_chip_fin = pchip((0:M-2)*dt,Average_FIN);
rho_chip_fin = @(t) ppval(pp_chip_fin,t);

%% 2.0 Reconstruction

% 2.1 Various ODEs
Adq = @(q) AdH(inv_se3(exp_se3(skew(q)))); % Convert right algebra element to left algebra element
dqdt_param = @(t,q) K(q,'inverse')*Adq(q)*rho_param(t); 
dqdt_fft = @(t,q) K(q,'inverse')*Adq(q)*rho_fit_FFT(t); 
dqdt_spline = @(t,q) K(q,'inverse')*Adq(q)*rho_fit_spline(t); 
dqdt_param_fin = @(t,q) K(q,'inverse')*Adq(q)*rho_param_Fin(t); 
dqdt_fft_fin = @(t,q) K(q,'inverse')*Adq(q)*rho_fit_FFT_Fin(t); 
dqdt_spline_fin = @(t,q) K(q,'inverse')*Adq(q)*rho_fit_spline_fin(t);
dqdt_chip = @(t,q) K(q,'inverse')*Adq(q)*rho_chip(t);
dqdt_chip_fin = @(t,q) K(q,'inverse')*Adq(q)*rho_chip_fin(t);

% 2.2 Solutions 
sol_param = ode45(@(t,q) dqdt_param(t,q), [0 T], zeros(6,1));
sol_fft = ode45(@(t,q) dqdt_fft(t,q), [0 T], zeros(6,1));
sol_spline = ode45(@(t,q) dqdt_spline(t,q), [0 T], zeros(6,1));
sol_param_fin = ode45(@(t,q) dqdt_param_fin(t,q), [0 T], zeros(6,1));
sol_fft_fin = ode45(@(t,q) dqdt_fft_fin(t,q), [0 T], zeros(6,1));
sol_spline_fin = ode45(@(t,q) dqdt_spline_fin(t,q), [0 T], zeros(6,1));
sol_chip = ode45(@(t,q) dqdt_chip(t,q), [0 T], zeros(6,1));
sol_chip_fin = ode45(@(t,q) dqdt_chip_fin(t,q), [0 T], zeros(6,1));

H_param = @(t) exp_se3(skew(deval(sol_param,t)));
H_fft = @(t) exp_se3(skew(deval(sol_fft,t)));
H_spline = @(t) exp_se3(skew(deval(sol_spline,t)));
H_param_fin = @(t) exp_se3(skew(deval(sol_param_fin,t)));
H_fft_fin = @(t) exp_se3(skew(deval(sol_fft_fin,t)));
H_spline_fin = @(t) exp_se3(skew(deval(sol_spline_fin,t)));
H_chip = @(t) exp_se3(skew(deval(sol_chip,t)));
H_chip_fin = @(t) exp_se3(skew(deval(sol_chip_fin,t)));

% 2.3 Reconstructed trajectories

l = 1; % chosen trajectory
k = l;
Sk_param = zeros(N*3, M);
Sk_param(:,1) = X0(:,k); % Initial condition
for j = 2:M
    Sk_param(:,j) = Phi(H_param(dt*(j-1)),Sk_param(:,1));
end

Sk_fft = zeros(N*3, M);
Sk_fft(:,1) = X0(:,k); % Initial condition
for j = 2:M
    Sk_fft(:,j) = Phi(H_fft(dt*(j-1)),Sk_fft(:,1));
end

Sk_spline = zeros(N*3, M);
Sk_spline(:,1) = X0(:,k); % Initial condition
for j = 2:M
    Sk_spline(:,j) = Phi(H_spline(dt*(j-1)),Sk_spline(:,1));
end

Sk_param_fin = zeros(N*3, M);
Sk_param_fin(:,1) = X0(:,k); % Initial condition
for j = 2:M
    Sk_param_fin(:,j) = Phi(H_param_fin(dt*(j-1)),Sk_param_fin(:,1));
end

Sk_fft_fin = zeros(N*3, M);
Sk_fft_fin(:,1) = X0(:,k); % Initial condition
for j = 2:M
    Sk_fft_fin(:,j) = Phi(H_fft_fin(dt*(j-1)),Sk_fft_fin(:,1));
end

Sk_spline_fin = zeros(N*3, M);
Sk_spline_fin(:,1) = X0(:,k); % Initial condition
for j = 2:M
    Sk_spline_fin(:,j) = Phi(H_spline_fin(dt*(j-1)),Sk_spline_fin(:,1));
end

Sk_chip = zeros(N*3, M);
Sk_chip(:,1) = X0(:,k); % Initial condition
for j = 2:M
    Sk_chip(:,j) = Phi(H_chip(dt*(j-1)),Sk_chip(:,1));
end

Sk_chip_fin = zeros(N*3, M);
Sk_chip_fin(:,1) = X0(:,k); % Initial condition
for j = 2:M
    Sk_chip_fin(:,j) = Phi(H_chip_fin(dt*(j-1)),Sk_chip_fin(:,1));
end

% 2.4 Errors (try arbitrary initial data not in testset?)
% 2.4.1 Over whole trajectory: 
Sk_GT = S_GT(:,(k-1)*M+1:k*M);
Ett_param = vecnorm(Sk_param - Sk_GT);
Ett_fft = vecnorm(Sk_fft - Sk_GT);
Ett_spline = vecnorm(Sk_spline - Sk_GT);
Ett_param_fin = vecnorm(Sk_param_fin - Sk_GT);
Ett_fft_fin = vecnorm(Sk_fft_fin - Sk_GT);
Ett_spline_fin = vecnorm(Sk_chip_fin - Sk_GT);
Ett_chip = vecnorm(Sk_chip - Sk_GT);
Ett_chip_fin = vecnorm(Sk_chip_fin - Sk_GT);


% 2.4.2 Between successive points:
Esp_param = zeros(M-2,1);
Esp_fft = zeros(M-2,1);
Esp_spline = zeros(M-2,1);
Esp_param_fin = zeros(M-2,1);
Esp_fft_fin = zeros(M-2,1);
Esp_spline_fin = zeros(M-2,1);
Esp_chip = zeros(M-2,1);
Esp_chip_fin = zeros(M-2,1);
for j = 2:M
    Esp_param(j-1) = vecnorm(Phi(exp_se3(dt*skew(rho_param((j-2)*dt))),Sk_GT(:,j-1)) - Sk_GT(:,j));
    Esp_fft(j-1) = vecnorm(Phi(exp_se3(dt*skew(rho_fit_FFT((j-2)*dt))),Sk_GT(:,j-1)) - Sk_GT(:,j));
    Esp_spline(j-1) = vecnorm(Phi(exp_se3(dt*skew(rho_fit_spline_smooth((j-2)*dt))),Sk_GT(:,j-1)) - Sk_GT(:,j));
    Esp_param_fin(j-1) = vecnorm(Phi(exp_se3(dt*skew(rho_param_Fin((j-2)*dt))),Sk_GT(:,j-1)) - Sk_GT(:,j));
    Esp_fft_fin(j-1) = vecnorm(Phi(exp_se3(dt*skew(rho_fit_FFT_Fin((j-2)*dt))),Sk_GT(:,j-1)) - Sk_GT(:,j));
    Esp_spline_fin(j-1) = vecnorm(Phi(exp_se3(dt*skew(rho_fit_spline_fin((j-2)*dt))),Sk_GT(:,j-1)) - Sk_GT(:,j));
    Esp_chip(j-1) = vecnorm(Phi(exp_se3(dt*skew(rho_chip((j-2)*dt))),Sk_GT(:,j-1)) - Sk_GT(:,j));
    Esp_chip_fin(j-1) = vecnorm(Phi(exp_se3(dt*skew(rho_chip_fin((j-2)*dt))),Sk_GT(:,j-1)) - Sk_GT(:,j));
end

%%
% 2.4.3 For multiple successive points
tmsp = tic;
Emsp_av = zeros(M-2,1);
Emsp = zeros(M-2,M-2);
for j = 2:M
    for jj = j:M
        Hjj = H_param((jj-1)*dt)/H_param((j-2)*dt);%H*exp_se3(dt*skew(rho_param((jj-2)*dt)));
        Emsp(j,jj-j+1) = vecnorm(Phi(Hjj,Sk_GT(:,j-1)) - Sk_GT(:,jj));
    end
end
Emsp_av = sum(Emsp,1);
for j = 1:M-1
    Emsp_av(j) = Emsp_av(j)/(M-j);
end
tmsp = toc(tmsp);


%% 3.0 Plotting
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


figure()
[U,SingVals,V] = svd(S,'econ');
plot(diag(SingVals),'.')
title('Singular values of pointcloud data')
xlabel('Index')
ylabel('Value')
legend(sprintf('Case for %i Particles, %i Trajectories',N,L))
grid on 
grid minor

figure()
[U1,SingVals1,V1] = svd(rho_Data,'econ');
plot(diag(SingVals1),'o')
title('Singular values of velocity-based reduced snapshots')
xlabel('Index')
ylabel('Value')
legend(sprintf('Case for %i Particles, %i Trajectories',N,L))
grid on 
grid minor

figure()
[U1,SingVals1,V1] = svd(rho_Data_Fin,'econ');
plot(diag(SingVals1),'o')
title('Singular values of velocity-free reduced snapshots')
xlabel('Index')
ylabel('Value')
legend(sprintf('Case for %i Particles, %i Trajectories',N,L))
grid on 
grid minor

figure()
T_tt = (0:M-1)*dt;
% plot(T_tt,Ett_param,T_tt,Ett_fft,T_tt,Ett_spline,T_tt,Ett_chip,T_tt,Ett_param_fin,T_tt,Ett_fft_fin,T_tt,Ett_spline_fin,T_tt,Ett_chip_fin,'LineWidth',2)
plot(T_tt,Ett_param,T_tt,Ett_spline,T_tt,Ett_chip,T_tt,Ett_param_fin,T_tt,Ett_spline_fin,T_tt,Ett_chip_fin,'LineWidth',2)

title('Error of MorLie: full trajectory')
xlabel('Time')
ylabel('Error')
% legend('Parametric Fit','Fourier Fit', 'Spline Fit','Hermite Fit', 'Parametric Fit: Finite twist', 'Fourier Fit: Finite twist', 'Spline Fit: Finite twist','Hermite Fit: Finite twist')
legend('Parametric Fit', 'Spline Fit','Hermite Fit', 'Parametric Fit: Finite twist', 'Spline Fit: Finite twist','Hermite Fit: Finite twist','Location','northwest')
grid on 
grid minor

figure()
T_sp = (0:M-2)*dt;
% plot(T_sp,Esp_param,T_sp,Esp_fft,T_sp,Esp_spline,T_sp,Esp_chip,T_sp,Esp_param_fin,T_sp,Esp_fft_fin,T_sp,Esp_spline_fin,T_sp,Esp_chip_fin,'LineWidth',2)
plot(T_sp,Esp_param,T_sp,Esp_spline,T_sp,Esp_chip,T_sp,Esp_param_fin,T_sp,Esp_spline_fin,T_sp,Esp_chip_fin,'LineWidth',2)
title('Error of MorLie: step ahead prediction')
xlabel('Time')
ylabel('Error')
% legend('Parametric Fit','Fourier Fit', 'Spline Fit','Hermite Fit', 'Parametric Fit: Finite twist', 'Fourier Fit: Finite twist', 'Spline Fit: Finite twist','Hermite Fit: Finite twist')
legend('Parametric Fit', 'Spline Fit','Hermite Fit', 'Parametric Fit: Finite twist', 'Spline Fit: Finite twist','Hermite Fit: Finite twist','Location','northwest')
grid on 
grid minor


figure()
T_sp = (0:M-2)*dt;
% plot(T_sp,Esp_param,T_sp,Esp_fft,T_sp,Esp_spline,T_sp,Esp_chip,T_sp,Esp_param_fin,T_sp,Esp_fft_fin,T_sp,Esp_spline_fin,T_sp,Esp_chip_fin,'LineWidth',2)
plot(T_sp,Emsp_av,'LineWidth',2)
title('Average Error: multiple step ahead prediction')
xlabel('Time')
ylabel('Error')
% legend('Parametric Fit','Fourier Fit', 'Spline Fit','Hermite Fit', 'Parametric Fit: Finite twist', 'Fourier Fit: Finite twist', 'Spline Fit: Finite twist','Hermite Fit: Finite twist')
legend('Parametric Fit','Location','northwest')
grid on 
grid minor

% Noisy video
% fig = figure();
% clear Frames;
% Frames(M) = struct('cdata', [], 'colormap', []);
% day = datetime('now','Format','y-M-d');
% time = datetime('now','Format','H-m-s');
% video_filename = sprintf('Rigid_point_cloud_noise_%s_at_%s',day,time);
% Sl = reshape(S(:,(l-1)*M+1:l*M),[3, N*M]);
% xm = min(Sl(1,:)); xp = max(Sl(1,:));
% ym = min(Sl(2,:)); yp = max(Sl(2,:));
% zm = min(Sl(3,:)); zp = max(Sl(3,:));
% pad = abs(xp-xm)*0.05;
% % title('Reconstruction of rigid point cloud')
% % legend('Full $\gamma(t)$','Reconstructed $\bar{\gamma} = \Phi\big(g(t),\gamma_0\big)$',interpreter='latex')
% % legend('')
% for j = 1:M 
%     Xj_GT = reshape(S(:,(l-1)*M+j),[3 N]);
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


%% Functions

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

function [pp] = coefs_to_pp(pp,coefs)
    pp.coefs = coefs;
end