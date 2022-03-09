clear;clc;close;
addpath(genpath("./"));
%%
Ns = 20;
dt = 0.025;
Nsim = ceil(6 / dt);
%%
param.Nx = 12; % dimension x
param.Nu = 6; % dimension of u
param.Nt = 20; % horizon length
param.dt = dt; % step time
param.xmax = inf(param.Nx, 1); % state constraints
param.xmin = -param.xmax;

Ib = [1,0,0;
    0,2,0;
    0,0,3];
M = eye(3) * 1;
I = [Ib, zeros(3);...
    zeros(3),  M];

param.umax = ones(param.Nu, 1) * 4000; % input constraints
param.umin = -param.umax;
param.I = I; % body inertial tensor
%% generate reference trajectory
q0_ref = [1;0;0;0];
p0_ref = [0;0;0];
w0_ref = [0;0;1] * 1;
v0_ref = [1;0;0.1] * 2;
x0_ref = [q0_ref;p0_ref];
xid_ref = [w0_ref;v0_ref];

X = eye(4);

X_ref = x0_ref';
xi_ref = xid_ref';
for i = 1:Nsim
    xid_ref_rt = xid_ref';
    % you can try some time verying twists:
    %   xid_ref_rt(1) = sin(i / 20) * 2;
    %   xid_ref_rt(5) = cos(sqrt(i)) * 1;
    %   xid_ref_rt(6) = 1; % sin(sqrt(i)) * 1;
    Xi = [skew(xid_ref_rt(1:3)), xid_ref_rt(4:6)';...
        [0,0,0,0]];
    X = X * expm(Xi * dt);
    X_ref = [X_ref; [rotm2quat(X(1:3,1:3)), X(1:3,4)']];
    xi_ref = [xi_ref; xid_ref_rt];
end
param.X_ref = X_ref;
param.xi_ref = xi_ref;
%%
Logger1 = struct('X',zeros(13, Nsim - 11),'Err',zeros(6, Nsim - 12),'U',zeros(6, Nsim - 12),'Err2', zeros(Nsim - 12, 4));
Logger2 = Logger1;
Logger3 = Logger1;

parfor k = 1:Ns
    q0 = rand(4,1) * 2 - 1;
    q0 = q0 / norm(q0);
    p0 = (rand(3,1) - 0.5) * 2;
    w0 = (rand(3,1) - 0.5) * 0;
    v0 = (rand(3,1) - 0.5) * 0;

    Logger1(k) = sim_mpc(q0, p0, w0, v0, dt, Nsim, 1, param);
    Logger2(k) = sim_mpc(q0, p0, w0, v0, dt, Nsim, 2, param);
    Logger3(k) = sim_mpc(q0, p0, w0, v0, dt, Nsim, 3, param);
    k
end
%%

axis_font_size = 12;
lw_1 = 0.1;

x_lim = [0, 5.5];
figure(1)
subplot(3,2,1)
for k = 1:Ns
    plot(0:dt:(length(Logger1(k).Err2)-1) * dt, Logger1(k).Err2(:,1),'LineWidth',lw_1)
    hold on
end
title("Proposed MPC - Orientation error", "Interpreter","latex")
box on
grid on
xlim(x_lim)
set(gca,'FontSize',axis_font_size)

subplot(3,2,3)
for k = 1:Ns
    plot(0:dt:(length(Logger2(k).Err2)-1) * dt, Logger2(k).Err2(:,1),'LineWidth',lw_1)
    hold on
end
title("VBL-MPC - Orientation error", "Interpreter","latex")
ylabel("$\|Log(R_d^{-1}R)\|$", "Interpreter","latex")
box on
grid on
xlim(x_lim)
set(gca,'FontSize',axis_font_size)

subplot(3,2,5)
for k = 1:Ns
    plot(0:dt:(length(Logger3(k).Err2)-1) * dt, Logger3(k).Err2(:,1),'LineWidth',lw_1)
    hold on
end
title("Simplified - Orientation error", "Interpreter","latex")
xlabel("Time (s)", "Interpreter","latex")
box on
grid on
xlim(x_lim)
set(gca,'FontSize',axis_font_size)

subplot(3,2,2)
for k = 1:Ns
    plot(0:dt:(length(Logger1(k).Err2)-1) * dt, vecnorm(Logger1(k).Err2(:,2:4),2,2),'LineWidth',lw_1)

    hold on
end
title("Proposed MPC - Position error", "Interpreter","latex")
box on
grid on
xlim(x_lim)
set(gca,'FontSize',axis_font_size)

subplot(3,2,4)
for k = 1:Ns
    plot(0:dt:(length(Logger2(k).Err2)-1) * dt, vecnorm(Logger2(k).Err2(:,2:4),2,2),'LineWidth',lw_1)
    hold on
end
title("VBL-MPC - Position error", "Interpreter","latex")
ylabel("$\|p-p_d\|$", "Interpreter","latex")
box on
grid on
xlim(x_lim)
set(gca,'FontSize',axis_font_size)

subplot(3,2,6)
for k = 1:Ns
    plot(0:dt:(length(Logger2(k).Err2)-1) * dt, vecnorm(Logger3(k).Err2(:,2:4),2,2),'LineWidth',lw_1)
    hold on
end
title("Simplified - Position error", "Interpreter","latex")
xlabel("Time (s)", "Interpreter","latex")
box on
grid on
xlim(x_lim)
set(gca,'FontSize',axis_font_size)

%%
figure(3)
eR_hist = zeros(Ns, 3);
eP_hist = zeros(Ns, 3);
for k = 1:Ns
    eR_hist(k, 1) = sum(Logger1(k).Err2(:,1));
    eR_hist(k, 2) = sum(Logger2(k).Err2(:,1));
    eR_hist(k, 3) = sum(Logger3(k).Err2(:,1));

    eP_hist(k, 1) = sum(vecnorm(Logger1(k).Err2(:,2:4)));
    eP_hist(k, 2) = sum(vecnorm(Logger2(k).Err2(:,2:4)));
    eP_hist(k, 3) = sum(vecnorm(Logger3(k).Err2(:,2:4)));
end
figure(3)
subplot(2,1,1)
hist(eR_hist(:,[1,2, 3]))
title("Histogram of orientation error")
subplot(2,1,2)
hist(eP_hist)
title("Histogram of position error")
legend({"Proposed MPC", "VBL-MPC", "Simplified"})
xlabel("Sum of absolute value of error")

% for k = 1:3
%     subplot(3, 1, k)
%     hist(Err_hist(:, k))
%     xlim([0. max(max(Err_hist))])
% end
%%
figure(4)
loggers = {Logger1, Logger2, Logger3};
for id = 1:3
    lw = 2;
    logger = loggers{id};
    if id == 1
        color = "r";
    end

    if id == 2
        color = "g";
    end

    if id == 3
        color = "b";
    end

    for k = 1
        plot3(logger(k).X(5, :), logger(k).X(6, :), logger(k).X(7, :), color+"-.", "LineWidth", lw)
        hold on
        %     plot3(logger(k).X(5, 1), logger(k).X(6, 1), logger(k).X(7, 1), "o","LineWidth",2)
        %     hold on
    end
end
lenn = 220;
plot3(X_ref(1:lenn,5), X_ref(1:lenn,6), X_ref(1:lenn,7),'-','color','k' ,"LineWidth",6)
daspect([1,1,1])
box on
grid on
legend({"Proposed MPC", "VBL-MPC", "Proposed MPc with simplified matrix", "Rererence"}, "interpreter", "latex")
xlabel("x", "interpreter", "latex")
ylabel('y', "interpreter", "latex")
zlabel('z', "interpreter", "latex")
%%
figure(5)
axis_font_size = 16;
k = 2;
for tt = 1:3
    subplot(1,3,tt)
    set(gca,'FontSize',18)
    if tt == 1
        logger = Logger1;
        title("Proposed MPC")
    end

    if tt == 2
        logger = Logger2;
        title("VBL MPC")
    end

    if tt == 3
        logger = Logger3;
        title("Simplified version")
    end


    h1 = plot3(X_ref(1:lenn,5), X_ref(1:lenn,6), X_ref(1:lenn,7),'-.','color','k' ,"LineWidth",2);
    hold on
    h2 = plot3(logger(k).X(5, :), logger(k).X(6, :), logger(k).X(7, :), "-",'color',[1,1,1]/2,"LineWidth", 3);

    arraw = 0.3;
    lw = 2;
    for kk = 1:30:length(logger(k).X)
        R = quat2rotm(logger(k).X(1:4,kk)');
        p = logger(k).X(5:7, kk);
        px = p + R(:,1) * arraw;
        py = p + R(:,2) * arraw;
        pz = p + R(:,3) * arraw;
        plot3([p(1), px(1)], [p(2), px(2)], [p(3), px(3)], "r-","LineWidth",lw)
        hold on
        plot3([p(1), py(1)], [p(2), py(2)], [p(3), py(3)], "g-","LineWidth",lw)
        hold on
        plot3([p(1), pz(1)], [p(2), pz(2)], [p(3), pz(3)], "b-","LineWidth",lw)
        hold on
    end
    box on
    grid on

    lw = 2;
    for kk = 1:30:length(X_ref) - 30
        R = quat2rotm(X_ref(kk,1:4));
        p = X_ref(kk,5:7)';
        px = p + R(:,1) * arraw;
        py = p + R(:,2) * arraw;
        pz = p + R(:,3) * arraw;
        plot3([p(1), px(1)], [p(2), px(2)], [p(3), px(3)], "-.","color",'r', "LineWidth",lw)
        hold on
        plot3([p(1), py(1)], [p(2), py(2)], [p(3), py(3)], "-.","color",'g',"LineWidth",lw)
        hold on
        plot3([p(1), pz(1)], [p(2), pz(2)], [p(3), pz(3)], "-.","color",'b',"LineWidth",lw)
        hold on
    end

    if tt == 1
        title("Proposed MPC")
    end

    if tt == 2
        title("VBL MPC")
    end

    if tt == 3
        title("Simplified version")
    end

    daspect([1,1,1])
    xlabel("$x$", "interpreter", "latex")
    ylabel('$y$', "interpreter", "latex")
    zlabel('$z$', "interpreter", "latex")
    % ylim([-0.2,4])
    set(gca,'FontSize',axis_font_size)
end
%%
legend([h1, h2], {"Reference (dashed)", "System Reponse (solid)"}, "interpreter", "latex")