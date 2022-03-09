function logger = sim_mpc(q0, p0, w0, v0, dt, Nsim, mode, param)
X_ref = param.X_ref;
xi_ref = param.xi_ref;

if mode ~= 2
    param.P = diag([1,1,1,...
        10,10,10,...
        1,1,1,...
        1,1,1]) * 10;
    param.Q = diag([1,1,1,...
        10,10,10,...
        1,1,1,...
        1,1,1]) * 1;
    param.R = eye(6) * 1e-5;
else
    param.P = diag([10,10,10,...
        1,1,1,...
        1,1,1,...
        1,1,1]) * 10;
    param.Q = diag([10,10,10,...
        1,1,1,...
        1,1,1,...
        1,1,1]) * 1;
    param.R = eye(6) * 1e-5;
end
%%
% log the sates
x0 = [q0;p0;w0;v0];
X = x0;
U = [];
Err = [];
Err2 = [];

logger = struct();

for i = 1:Nsim - param.Nt
    X_ref_rt = X_ref(i,:);
    xi_ref_rt = xi_ref(i:i+param.Nt-1, :);
    if mode == 1
        u = eMPC(x0, X_ref_rt,  xi_ref_rt, param);
    elseif mode == 2
        u = gMPC(x0, X_ref_rt,  xi_ref_rt, param);
    elseif mode == 3
        u = eMPC_simplified(x0, X_ref_rt,  xi_ref_rt, param);
    else
        assert(1)
    end

    [t,y] = ode45(@(t,x)SE3Dyn(t,x,u,param.I), [0:param.dt/20:param.dt], x0); %  + 0*abs(randn) * param.dt / 2
    x0 = y(end,:)';


    Xd = [quat2rotm(X_ref_rt(1:4)), X_ref_rt(5:7)';...
        0,0,0,1];
    X0 = [quat2rotm(x0(1:4)'), x0(5:7);...
        0,0,0,1];
    Xerr = logm(X0^-1 * Xd);
    XXerr = X0^-1 * Xd;
    x00 = x0(2:7);
    x00(1:3) = [Xerr(3,2); Xerr(1,3); Xerr(2,1)]; %% R error
    x00(4:6) = Xerr(1:3,4); %% p error

    X = [X, x0];
    U = [U, u];
    Err = [Err, x00];
    Err2 = [Err2; [sqrt(sum(sum(logm(XXerr(1:3,1:3)).^2))/2), XXerr(1:3,4)'] ];
end
logger.X = X;
logger.U = U;
logger.Err = Err;
logger.Err2 = Err2;
end