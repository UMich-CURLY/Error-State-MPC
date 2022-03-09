function u = eMPC_casadi(x0, xd, xid, param)
Xd = [quat2rotm(xd(1:4)), xd(5:7)';...
    0,0,0,1];
X00 = [quat2rotm(x0(1:4)'), x0(5:7);...
    0,0,0,1];

% The inverse of Left error defined in paper is equivalent.
X00 =  X00^(-1) * Xd;
X0 = logm(X00);

p0 = [X0(3,2); X0(1,3); X0(2,1); X0(1:3,4);x0(8:end)]; %% initial error
%%
I = param.I;
dt = param.dt;
Nx = param.Nx;
%% use Casadi:
addpath('D:\MatlabPackage\casadi-windows-matlabR2016a-v3.4.5')
import casadi.*
%%
x = SX.sym('x', [param.Nx, param.Nt + 1]);
u = SX.sym('u', [param.Nu, param.Nt]);

g_dynamics = [];
g_dynamics_min = [];
g_dynamics_max = [];
J = 0;
for k = 1:param.Nt
    % xi_bar = I * xid(k, :)';
    xi_bar = I * x0(8:13);
    G = zeros(6,6);
    G(1:3,1:3) = skew(xi_bar(1:3));
    G(1:3,4:6) = skew(xi_bar(4:6));
    G(4:6,1:3) = skew(xi_bar(4:6));
    
    H = I^-1 *(coadjoint(x0(8:13)) * I + G);
    
    Ac = [-adjoint_(xid(k, :)),   -eye(6);...
        zeros(6), H];
    Bc = [zeros(6);...
        I^-1];
    
    b = -I^(-1) * G * x0(8:13);
    hc = [xid(k,:)'; b];
    
    Ad = expm(Ac * dt);
    Bd = (eye(param.Nx) * dt + Ac * dt^2 / 2 + Ac^2 * dt^3 / 6) * Bc;
    hd = (eye(param.Nx) * dt + Ac * dt^2 / 2 + Ac^2 * dt^3 / 6) * hc;
    
    X_k = x(:, k);
    X_kp1 = x(:, k+1);
    u_k = u(:, k);
    
    g_dynamics = [g_dynamics;...
        X_kp1 - Ad * X_k - Bd * u_k - hd];
    g_dynamics_min = [g_dynamics_min;...
        zeros(Nx,1)];
    g_dynamics_max = [g_dynamics_max;...
        zeros(Nx,1)];
    
    C = eye(12);
    C(7:12, 1:6) = -adjoint_(xid(k,:));
    d = zeros(12,1);
    d(7:12) = xid(k,:)';
    J = J + 0.5 * (C * x(:,k) - d)' * param.Q * (C * x(:,k) - d) + 0.5 * u(:,k)' * param.R * u(:,k);
end

P = idare(Ad, Bd, param.Q,param.R,[], []);
J = J + 0.5 * (C * X_kp1 - d)' * P * (C * X_kp1 - d);

g_boundary = [];
g_boundary = [g_boundary; x(:,1) - p0]; % initial error

g = [g_dynamics; g_boundary];
g_min = [g_dynamics_min; zeros(param.Nx, 1)];
g_max = [g_dynamics_max; zeros(param.Nx, 1)];

b_input_min = repmat(param.umin, [param.Nt, 1]);
b_input_max = repmat(param.umax, [param.Nt, 1]);

x_max = [b_input_max; inf(param.Nx * (param.Nt+1), 1)];
x_min = [b_input_min; -inf(param.Nx * (param.Nt+1), 1)];

input = [u(:); x(:)];
%%
nlp = struct('x',input, 'f',J, 'g', g);
% help nlpsol
S = nlpsol('S', 'ipopt', nlp);

x0 = randn(size(input));
r = S('x0',x0,'lbg',g_min,'ubg',g_max,'lbx',x_min,'ubx',x_max);
sol = full(r.x);
u = sol(1:6);
end
