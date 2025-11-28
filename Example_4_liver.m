
%% 0. Setup
addpath(genpath('./Auxilliary'));

NAmb = 3;
L = 1; % Number of trajectories
T = 116; % Time
Sred = csvread("../liver_tracking.csv");
M = size(Sred,2);
N = size(Sred,1)/2;
dt = T/M; % Time step

S = zeros(NAmb*N,M);
for j = 1:size(Sred,1)
    i = mod(j-1,2)+1;
    I = floor((j-1)/2)*3;
    S(I+i,:) = Sred(j,:);
end

X0 = S(:,1);
S_GT = S;

%% 1.0 MorLie
% 1.1 Pre-processing:
% 1.1.1 Fitting Infinitessimal Generators to Vector fields

S_red = zeros(N*NAmb,(M-1)*L); % remove final point 
Sp = zeros(N*NAmb,(M-1)*L); % difference between successive points (on the same trajectory)
for k = 1:L
    Sk = S(:,(k-1)*M+1:k*M);
    Sk_red = Sk(:,1:end-1);
    Spk = Sk(:,2:end) - Sk(:,1:end-1);
    S_red(:,(k-1)*(M-1)+1:k*(M-1)) = Sk_red;
    Sp(:,(k-1)*(M-1)+1:k*(M-1)) = Spk;
end


Time_FitInfGen = tic;
dim_G = 6;
rho_Data = zeros(dim_G,(M-1)*L);
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
rho_Data_Fin = zeros(dim_G,(M-1)*L);
g_Data_Fin = cell(L,M);
q_init = zeros(dim_G,1);
options = optimoptions('fsolve','Display','none','MaxIterations',10,'Algorithm','levenberg-marquardt');
for k = 1:L
    g_Data_Fin{k,1} = eye(4);
    for j = 1:(M-1)
        ind = (k-1)*(M-1)+j;
        X = S_red(:,ind);
        dX = Sp(:,ind);
        X_next = X+dX; 
        rho_Data_Fin(:,ind) = unskew(log_se3(LM2_SE3(@(G,X) Phi(G,X),eye(4),X,X_next,2,0))); %fsolve(@(q)dist_Phi(exp_se3(skew(q)),X,X_next),q_init,options); %fminsearch(@(q)dist_Phi(exp_se3(skew(q)), X, X_next),q0); 
        q_init = rho_Data_Fin(:,ind);
        g_Data_Fin{k,j+1} = exp_se3(skew(rho_Data_Fin(:,ind)))*g_Data_Fin{k,j};
    end
end
rho_Data_Fin = rho_Data_Fin/dt;
Time_FitFinGen = toc(Time_FitFinGen);

%% 1.2 Fitting reduced vector field to rho_Data

T_red = repmat((0:M-2)*dt,1,L);

% 1.2.2 Truncated Fourier series to rho_Data
Average = zeros(dim_G,M-1);
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

% 1.2.6 Truncated Fourier series to rho_Data_Fin
Average_FIN = zeros(dim_G,M-1);
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

%% 2.0 Reconstruction: need to code AdG

% 2.1 Various ODEs
Adq = @(q) AdH(inv_se3(exp_se3(skew(q)))); % Convert right algebra element to left algebra element
dqdt_fft = @(t,q) K(q,'inverse')*Adq(q)*rho_fit_FFT(t); 
dqdt_spline = @(t,q) K(q,'inverse')*Adq(q)*rho_fit_spline(t); 
dqdt_fft_fin = @(t,q) K(q,'inverse')*Adq(q)*rho_fit_FFT_Fin(t); 
dqdt_spline_fin = @(t,q) K(q,'inverse')*Adq(q)*rho_fit_spline_fin(t);
dqdt_chip = @(t,q) K(q,'inverse')*Adq(q)*rho_chip(t);
dqdt_chip_fin = @(t,q) K(q,'inverse')*Adq(q)*rho_chip_fin(t);

% 2.2 Solutions 
sol_fft = ode45(@(t,q) dqdt_fft(t,q), [0 T], zeros(dim_G,1));
sol_spline = ode45(@(t,q) dqdt_spline(t,q), [0 T], zeros(dim_G,1));
sol_fft_fin = ode45(@(t,q) dqdt_fft_fin(t,q), [0 T], zeros(dim_G,1));
sol_spline_fin = ode45(@(t,q) dqdt_spline_fin(t,q), [0 T], zeros(dim_G,1));
sol_chip = ode45(@(t,q) dqdt_chip(t,q), [0 T], zeros(dim_G,1));
sol_chip_fin = ode45(@(t,q) dqdt_chip_fin(t,q), [0 T], zeros(dim_G,1));


H_fft = @(t) exp_se3(skew(deval(sol_fft,t)));
H_spline = @(t) exp_se3(skew(deval(sol_spline,t)));
H_fft_fin = @(t) exp_se3(skew(deval(sol_fft_fin,t)));
H_spline_fin = @(t) exp_se3(skew(deval(sol_spline_fin,t)));
H_chip = @(t) exp_se3(skew(deval(sol_chip,t)));
H_chip_fin = @(t) exp_se3(skew(deval(sol_chip_fin,t)));

% 2.3 Reconstructed trajectories

l = 1; % chosen trajectory
k = l;

Sk_fft = zeros(N*NAmb, M);
Sk_fft(:,1) = X0(:,k); % Initial condition
for j = 2:M
    Sk_fft(:,j) = Phi(H_fft(dt*(j-1)),Sk_fft(:,1));
end

Sk_spline = zeros(N*NAmb, M);
Sk_spline(:,1) = X0(:,k); % Initial condition
for j = 2:M
    Sk_spline(:,j) = Phi(H_spline(dt*(j-1)),Sk_spline(:,1));
end

Sk_fft_fin = zeros(N*NAmb, M);
Sk_fft_fin(:,1) = X0(:,k); % Initial condition
for j = 2:M
    Sk_fft_fin(:,j) = Phi(H_fft_fin(dt*(j-1)),Sk_fft_fin(:,1));
end

Sk_spline_fin = zeros(N*NAmb, M);
Sk_spline_fin(:,1) = X0(:,k); % Initial condition
for j = 2:M
    Sk_spline_fin(:,j) = Phi(H_spline_fin(dt*(j-1)),Sk_spline_fin(:,1));
end

Sk_chip = zeros(N*NAmb, M);
Sk_chip(:,1) = X0(:,k); % Initial condition
for j = 2:M
    Sk_chip(:,j) = Phi(H_chip(dt*(j-1)),Sk_chip(:,1));
end

Sk_chip_fin = zeros(N*NAmb, M);
Sk_chip_fin(:,1) = X0(:,k); % Initial condition
for j = 2:M
    Sk_chip_fin(:,j) = Phi(H_chip_fin(dt*(j-1)),Sk_chip_fin(:,1));
end

S_best = zeros(N*NAmb, M);
Sk_best(:,1) = X0(:,k); % Initial condition
for j = 2:M
    Sk_best(:,j) = Phi(g_Data_Fin{k,j},Sk_best(:,1));
end


% 2.4 Errors (try arbitrary initial data not in testset?)
% 2.4.1 Over whole trajectory: 
Sk_GT = S;
Ett_fft = vecnorm(Sk_fft - Sk_GT);
Ett_spline = vecnorm(Sk_spline - Sk_GT);
Ett_fft_fin = vecnorm(Sk_fft_fin - Sk_GT);
Ett_spline_fin = vecnorm(Sk_chip_fin - Sk_GT);
Ett_chip = vecnorm(Sk_chip - Sk_GT);
Ett_chip_fin = vecnorm(Sk_chip_fin - Sk_GT);
Ett_best = vecnorm(Sk_best - Sk_GT);

Ett_fft_fin_inf = vecnorm(Sk_fft_fin - Sk_GT,inf);
Ett_spline_fin_inf = vecnorm(Sk_chip_fin - Sk_GT,inf);
Ett_chip_fin_inf = vecnorm(Sk_chip_fin - Sk_GT,inf);
Ett_best_inf = vecnorm(Sk_best - Sk_GT,inf);


% 2.4.2 Between successive points:
Esp_fft = zeros(M-2,1);
Esp_spline = zeros(M-2,1);
Esp_fft_fin = zeros(M-2,1);
Esp_spline_fin = zeros(M-2,1);
Esp_chip = zeros(M-2,1);
Esp_chip_fin = zeros(M-2,1);

Esp_spline_fin_inf = zeros(M-2,1);
Esp_chip_fin_inf = zeros(M-2,1);
Esp_best = zeros(M-2,1);
Esp_best_inf = zeros(M-2,1);


for j = 2:M
    Esp_fft(j-1) = vecnorm(Phi(exp_se3(dt*skew(rho_fit_FFT((j-2)*dt))),Sk_GT(:,j-1)) - Sk_GT(:,j));
    Esp_spline(j-1) = vecnorm(Phi(exp_se3(dt*skew(rho_fit_spline_smooth((j-2)*dt))),Sk_GT(:,j-1)) - Sk_GT(:,j));
    Esp_fft_fin(j-1) = vecnorm(Phi(exp_se3(dt*skew(rho_fit_FFT_Fin((j-2)*dt))),Sk_GT(:,j-1)) - Sk_GT(:,j));
    Esp_spline_fin(j-1) = vecnorm(Phi(exp_se3(dt*skew(rho_fit_spline_fin((j-2)*dt))),Sk_GT(:,j-1)) - Sk_GT(:,j));
    Esp_chip(j-1) = vecnorm(Phi(exp_se3(dt*skew(rho_chip((j-2)*dt))),Sk_GT(:,j-1)) - Sk_GT(:,j));
    Esp_chip_fin(j-1) = vecnorm(Phi(exp_se3(dt*skew(rho_chip_fin((j-2)*dt))),Sk_GT(:,j-1)) - Sk_GT(:,j));
    Esp_spline_fin_inf(j-1) = vecnorm(Phi(exp_se3(dt*skew(rho_fit_spline_fin((j-2)*dt))),Sk_GT(:,j-1)) - Sk_GT(:,j),inf);
    Esp_chip_fin_inf(j-1) = vecnorm(Phi(exp_se3(dt*skew(rho_chip_fin((j-2)*dt))),Sk_GT(:,j-1)) - Sk_GT(:,j),inf);
    Esp_best(j-1) = vecnorm(Phi(g_Data_Fin{k,j}/g_Data_Fin{k,j-1},Sk_GT(:,j-1)) - Sk_GT(:,j));
    Esp_best_inf(j-1) = vecnorm(Phi(g_Data_Fin{k,j}/g_Data_Fin{k,j-1},Sk_GT(:,j-1)) - Sk_GT(:,j),inf);
end

%%
% 2.4.3 For multiple successive points
% Emsp_av = zeros(M-2,1);
% Emsp = zeros(M-2,M-2);
% for j = 2:M
%     for jj = j:M
%         Hjj = H_chip((jj-1)*dt)/H_chip((j-2)*dt);%H*exp_se3(dt*skew(rho_param((jj-2)*dt)));
%         Emsp(j,jj-j+1) = vecnorm(Phi(Hjj,Sk_GT(:,j-1)) - Sk_GT(:,jj));
%     end
% end
% Emsp_av = sum(Emsp,1);
% for j = 1:M-1
%     Emsp_av(j) = Emsp_av(j)/(M-j);
% end

%% 3.0 Plotting

day = datetime('now','Format','y-M-d');
time = datetime('now','Format','H-m-s');
FolderName = sprintf('Liver_%s_at_%s',day,time);
mkdir(FolderName);

fig = figure();
clear Frames;
skip = 1;
Frames(M/skip) = struct('cdata', [], 'colormap', []);
day = datetime('now','Format','y-M-d');
time = datetime('now','Format','H-m-s');
% video_filename = sprintf('liver_tracking_%s_at_%s',day,time);
Sl = reshape(S(:,(l-1)*M+1:l*M),[3, N*M]);
xm = min(Sl(1,:)); xp = max(Sl(1,:));
ym = min(Sl(2,:)); yp = max(Sl(2,:));
zm = min(Sl(3,:)); zp = max(Sl(3,:));
pad = abs(xp-xm)*0.05;
% % title('Reconstruction of rigid point cloud')
% % legend('Full $\gamma(t)$','Reconstructed $\bar{\gamma} = \Phi\big(g(t),\gamma_0\big)$',interpreter='latex')
% % legend('')
% ind = 0;
% for j = 1:skip:M
%     ind = ind+1;
%     Xj_GT = reshape(S_GT(:,(l-1)*M+j),[3 N]);
%     Xj_fit = reshape(Sk_best(:,j),[3 N]);
%     x_GT = Xj_GT(1,:); y_GT = Xj_GT(2,:); z_GT = Xj_GT(3,:);
%     x_fit = Xj_fit(1,:); y_fit = Xj_fit(2,:); z_fit = Xj_fit(3,:);
%     scatter3(x_GT,y_GT,z_GT,'o');
%     hold on
%     scatter3(x_fit,y_fit,z_fit,'*');
%     hold off
%     view(0,-90);
%     title('Reconstruction of rigid point cloud')
%     legend('Full $\gamma(t)$','Reconstructed $\bar{\gamma} = \Phi\big(g(t),\gamma_0\big)$',interpreter='latex',location='northeast')
%     xlim([xm-pad xp+pad]); ylim([ym-pad yp+pad]); zlim([zm-pad zp+pad]);
%     Frames(ind) = getframe(fig);
% end
% video = VideoWriter(sprintf('%s.avi',video_filename));
% video.FrameRate = 1/dt/skip;
% open(video)
% writeVideo(video,Frames(2:end));
% close(video)

ind = 0;
for j = round([1,M/4,M/2,M/4*3,M])
    fig = figure();
    ind = ind+1;
    name = sprintf('ReconstructionSnapshot_%iL.svg',ind);
    name = strcat(FolderName,'/',name);
    Xj_GT = reshape(S_GT(:,(l-1)*M+j),[3 N]);
    Xj_fit = reshape(Sk_best(:,j),[3 N]);
    x_GT = Xj_GT(1,:); y_GT = Xj_GT(2,:); z_GT = Xj_GT(3,:);
    x_fit = Xj_fit(1,:); y_fit = Xj_fit(2,:); z_fit = Xj_fit(3,:);
    scatter3(x_GT,y_GT,z_GT,'o');
    hold on
    scatter3(x_fit,y_fit,z_fit,'*');
    hold off
    view(0,-90);
    title(sprintf('T = %.2f s', j*dt))
    legend('Full $P_i(t)$','Reconstructed $\bar{P}_i(t) = \Phi\big(g(t), P_{i,0}\big)$',interpreter='latex',location='southeast')
    xlabel('mm');
    ylabel('mm');
    fontsize(17,"points");
    xlim([xm-pad xp+pad]); ylim([ym-pad yp+pad]); zlim([zm-pad zp+pad]);
    saveas(fig,name);
end

name = sprintf('SVD.fig');
name = strcat(FolderName,'/',name);
name2 = sprintf('SVD.svg');
name2 = strcat(FolderName,'/',name2);
[U,SingVals,V] = svd(S,'econ');
semilogy(diag(SingVals),'.')
title('Singular values of $S$')
xlabel('Index')
ylabel('Value')
legend('Singular values');%printf('Case for %i Particles, %i Trajectories',N,L))
fontsize(17,"points");
grid on 
grid minor
saveas(fig,name);
saveas(fig,name2);

name = sprintf('SVD_VF.fig');
name = strcat(FolderName,'/',name);
name2 = sprintf('SVD_VF.svg');
name2 = strcat(FolderName,'/',name2);
fig = figure();
[U1,SingVals1,V1] = svd(rho_Data_Fin,'econ');
semilogy(diag(SingVals1),'o')
title('Singular values of $S_{se(3)}$')
xlabel('Index')
ylabel('Value')
legend('Velocity-free')
fontsize(17,"points");
grid on 
grid minor
saveas(fig,name);
saveas(fig,name2);

fig  = figure();
T_tt = (0:M-1)*dt;
plot(T_tt,Ett_spline,T_tt,Ett_chip,T_tt,Ett_spline_fin,T_tt,Ett_chip_fin,'LineWidth',2)


title('Error of MorLie: full trajectory')
xlabel('Time')
ylabel('Error')
legend( 'Spline Fit','Hermite Fit', 'Spline Fit: Finite twist','Hermite Fit: Finite twist','Location','northwest')
fontsize(17,"points");
grid on 
grid minor
name = sprintf('TT_Error.svg');
name = strcat(FolderName,'/',name);
saveas(fig,name);


fig = figure();
T_sp = (0:M-2)*dt;
plot(T_sp,Esp_spline,T_sp,Esp_chip,T_sp,Esp_spline_fin,T_sp,Esp_chip_fin,'LineWidth',2)
title('Step ahead error')
xlabel('Time')
ylabel('Error')
legend( 'Spline Fit','Hermite Fit','Spline Fit: Finite twist','Hermite Fit: Finite twist','Location','northwest')
fontsize(17,"points");
grid on 
grid minor
name = sprintf('SP_Error.svg');
name = strcat(FolderName,'/',name);
saveas(fig,name);

fig = figure();
plot(T_tt,Ett_spline_fin/N,T_tt,Ett_chip_fin/N,T_tt,Ett_spline_fin_inf,T_tt,Ett_chip_fin_inf,'LineWidth',2)

title('Full trajectory error')
xlabel('Time')
ylabel('Error')
legend( 'Spline Fit: average error','Hermite Fit: average error','Spline Fit: max error','Hermite Fit: max error','Location','northwest')
fontsize(17,"points");
grid on 
grid minor
name = sprintf('SP_Error_InfNorm.svg');
name = strcat(FolderName,'/',name);
saveas(fig,name);

fig = figure();
plot(T_sp,Esp_spline_fin/N,T_sp,Esp_chip_fin/N,T_sp,Esp_spline_fin_inf,T_sp,Esp_chip_fin_inf,'LineWidth',2)

title('Step ahead error')
xlabel('Time')
ylabel('Error')
legend( 'Spline Fit: average error','Hermite Fit: average error','Spline Fit: max error','Hermite Fit: max error','Location','northwest')
fontsize(17,"points");
grid on 
grid minor
name = sprintf('TT_Error_InfNorm.svg');
name = strcat(FolderName,'/',name);
saveas(fig,name);

fig = figure();
plot(T_tt,Ett_best/N,T_tt,Ett_best_inf,'LineWidth',2)

title('Full trajectory error')
xlabel('Time in s')
ylabel('Error in mm')
legend( 'Average error','Max error','Location','northwest')
fontsize(17,"points");
grid on 
grid minor
name = sprintf('TT_Error_Best.svg');
name = strcat(FolderName,'/',name);
saveas(fig,name);


fig = figure();
plot(T_sp,Esp_best/N,T_sp,Esp_best_inf,'LineWidth',2)

title('Step ahead error')
xlabel('Time in  s')
ylabel('Error in mm')
legend( 'Average error','Max error','Location','northwest')
fontsize(17,"points");
grid on 
grid minor
name = sprintf('SP_Error_Best.svg');
name = strcat(FolderName,'/',name);
saveas(fig,name);

% figure()
% T_sp = (0:M-2)*dt;
% % plot(T_sp,Esp_param,T_sp,Esp_fft,T_sp,Esp_spline,T_sp,Esp_chip,T_sp,Esp_param_fin,T_sp,Esp_fft_fin,T_sp,Esp_spline_fin,T_sp,Esp_chip_fin,'LineWidth',2)
% plot(T_sp,Emsp_av,'LineWidth',2)
% title('Average Error: multiple step ahead prediction')
% xlabel('Time')
% ylabel('Error')
% % legend('Parametric Fit','Fourier Fit', 'Spline Fit','Hermite Fit', 'Parametric Fit: Finite twist', 'Fourier Fit: Finite twist', 'Spline Fit: Finite twist','Hermite Fit: Finite twist')
% legend('Parametric Fit','Location','northwest')
% grid on 
% grid minor

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

% Jacobian mapping components of A to X_A (dimension 3*N by dim_G)
function [PhiJac] = JacPhi(X)
    [N3,~] = size(X); N = N3/3;
    Ai = eye(6);
    PhiJac = zeros(N*3,6);
    for i = 1:6
        X_A = skew(Ai(:,i))*[reshape(X,[3 N]); ones(1,N)];
        PhiJac(:,i) = reshape(X_A(1:3,1:N),[N*3 1]);
    end
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


function [error] = dist_Phi(H,X,Xp)
    error = norm(Xp - Phi(H,X));
end

function [error] = dist_Phi_vec(H,X,Xp)
    error = (Xp - Phi(H,X)).^2;
end

function [pp] = coefs_to_pp(pp,coefs)
    pp.coefs = coefs;
end

function [Atilde] = skewgl(A)
    Atilde = [reshape(A,[3, 4]); [0, 0, 0, 0]];
end

function [A] = unskewgl(Atilde)
    A = reshape(Atilde(1:3,:),[12,1]);
end

function [G] = expgl(A) % Bad exp, only locally valid (small A)
    G = eye(4) + skewgl(A);
end

function [A] = loggl(G) % Bad log, only locally valid (small A)
    A = unskewgl(G - eye(4));
end