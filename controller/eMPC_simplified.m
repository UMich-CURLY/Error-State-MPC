function u = eMPC_simplified(x0, xd, xid, param)
% [φ, θ, ψ, x, y, z, ω1, ω2, ω3, vx, vy, vz]
R = quat2rotm(x0(1:4)');
I = param.I;
dt = param.dt;
Xd = [quat2rotm(xd(1:4)), xd(5:7)';...
      0,0,0,1];
X00 = [quat2rotm(x0(1:4)'), x0(5:7);...
      0,0,0,1];

X00 =  X00^(-1) * Xd; %% The error
% X00(1:3,end) = 0;
X0 = logm(X00);

Ad = [eye(6), -eye(6) * dt;...
      zeros(6), eye(6)];
Bd = [-1/2 * I^(-1) * dt^2;...
      I^(-1) * dt];

b = [xid(1,:)' * dt;...
     I^(-1) * (coadjoint(x0(8:13)) * I * x0(8:end)) * dt];
x0 = [X0(3,2); X0(1,3); X0(2,1); X0(1:3,4);x0(8:end)]; %% error

param.Ad = eye(6);

param.P = idare(Ad, Bd,param.Q,param.R,[], []);
[M, q] = eMPCCost_simplified(param.Q, param.R, param.P, xid, param);

[A, bmin, bmax] = MPCConstraints(Ad, Bd, b, x0, param);
solver = osqp;
solver.setup(M, q, A, bmin, bmax,'verbose', 0);
sol = solver.solve;
u = sol.x((param.Nt+1)*param.Nx+1:(param.Nt+1)*param.Nx+param.Nu);
end
% 
% I = param.I;
% dt = param.dt;
% 
% Ad = [eye(6), -eye(6) * dt;...
%      zeros(6), eye(6)];
% Bd = [-1/2 * I^(-1) * dt^2;...
%      I^(-1) * dt];
% b = I^(-1) * coadjoint(x0) * x0(7:end);
% b = [-b * dt^2 / 2 + AdjointMat(x0) * I * xid * dt;b * dt];
% 
% M = MPCCost(param.Q, param.R, param.P, param);
% X0 = [eul2rotm([x0(1), x0(2), x0(3)]), x0(4:6);...
%       0,0,0,1];
% X0 = logm(X0);
% x0(1:6) = [X0(3,2); X0(1,3); X0(2,1); X0(1:3,4)]; %% error
% 
% [A, bmin, bmax] = MPCConstraints(Ad, Bd, b, x0, param);
% solver = osqp;
% solver.setup(M, zeros(size(M,1),1), A, bmin, bmax);
% sol = solver.solve;
% u = sol.x((param.Nt+1)*param.Nx+1:(param.Nt+1)*param.Nx+param.Nu);
% end