%% 1. Setup

%% 1.1 Full order system (FOM)
a = 3; b = 5*2*3; omega = 1;                                                % large oscillations
%a = 1; b = 50*2*pi;                                                        % small oscillations

dxdt = @(x) [a*cos(b*x(2));omega];                                          % dynamics of radial oscillator
T = 2*pi;                                                                   % upper time for simulation

%% 1.2 Reduced order system (ROM)
Phi = @(g,x) [x(1)+g*0;x(2)+g];                                             % G = SO(2) in representation as real number, elements act on (r,theta) by adding to angle
rho = @(x) omega; %[0, -1; 1, 0];                                               
id_G = 0; 
%exp_G = @(A) [cos(A), -sin(A); sin(A), cos(A)];


%% 2. Solving ODEs:
options = odeset('AbsTol',1e-10,'RelTol',1e-10,'MaxStep',0.1);

n = 3;                                                                      % number of ICs
lower = 0.5;                                                                % step for ICs
x0s = [lower:lower:n*lower; zeros(1,n)];

xn = cell(n,1);                                                             % full order states
gn = cell(n,1);                                                             % reduced order states
x_barn = cell(n,1);                                                         % approximate states

for i = 1:n
    %% 2.1 Solving FOM
    %fprintf("Simulating FOM ... \n")
    tic
    x0 = x0s(:,i);
    x = ode45(@(t,x)dxdt(x),[0 T], x0, options);
    xn{i} = x.y;
    t_FOM(i) = toc;

    %% 2.2 Solving ROM
    %fprintf("Simulating ROM ... \n")
    tic
    g = ode45(@(t,g)rho(Phi(g,x0)),[0 T], id_G, options);
    gn{i} = g.y;
    x_barn{i} = Phi(gn{i},x0);
    t_ROM(i) = toc;
end

%% 3. Plotting
polar_to_cart = @(x) [x(1,:).*cos(x(2,:));x(1,:).*sin(x(2,:))];

figure()
hold on
grid on
ticks = -1.6:0.4:1.6;

for i = 1:n
    cart_FOM = polar_to_cart(xn{i});
    cart_ROM = polar_to_cart(x_barn{i});  
    plot(cart_FOM(1,:),cart_FOM(2,:), 'linewidth',2, 'Color','#A2142F')
    plot(cart_ROM(1,:),cart_ROM(2,:), 'linewidth',2, 'Color','#77AC30')
    %legendstr = legendstr + ',' + 'FOM','ROM'
end

hold off
title('Radial oscillator',Interpreter ='latex')
xlabel('$x$-Position',Interpreter='latex')
ylabel('$y$-Position',Interpreter='latex')
legend('FOM','ROM',Interpreter='latex')

xticks(ticks)
xlim([-1.6, 1.6])
yticks(ticks)
ylim([-1.6, 1.6])