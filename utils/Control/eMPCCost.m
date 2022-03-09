function [M, q] = eMPCCost(Q, R, P, xid, param)
% min x'Qx + (xi-xid)'Q(xi-xid) + u'Ru
% xkp1 = Ak xk + Bk uk + bk, k = 0,1, ..., N-1
% x0 = x_init,
% umin < uk < umax, k = 0, 1, ..., N-1
M = Q * 0;
q = zeros((param.Nt+1)*param.Nx + param.Nt * param.Nu,1);
for k = 1:param.Nt-1
    %%% trivial one
    %     M = blkdiag(M, Q);
    %     w = diag(Q);
    %     q((k-1)*param.Nx+7:(k)*param.Nx) = - w(7:end) .* (param.Ad * xid(k,:)'); %
    %% consider the transport map
    C = eye(12);
    C(7:12, 1:6) = -adjoint_(xid(k,:));
    M = blkdiag(M, C' * Q * C);
    b = zeros(12,1);
    b(7:12) = xid(k,:)';
    b = C' * Q * b;
    q((k-1)*param.Nx+1:(k)*param.Nx) = -b;
end

k = param.Nt;
C = eye(12);
C(7:12, 1:6) = -adjoint_(xid(k,:));
b = zeros(12,1);
b(7:12) = xid(k,:)';
b = C' * P * b;
M = blkdiag(M, C'*P*C);
q((k-1)*param.Nx+1:(k)*param.Nx) = -b; % param.Ad *
for k = 1:param.Nt
    M = blkdiag(M, R);
end
end