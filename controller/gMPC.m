function u = gMPC(x0, xd, xid, param)
    R = quat2rotm(x0(1:4)');
    I = param.I;
    dt = param.dt;
    
    Rd = quat2rotm(xd(1:4));
    pd = xd(5:7)';
    wd = xid(1, 1:3)';
    vd = xid(1, 4:6)';
    
    R0 = quat2rotm(x0(1:4)');
    p0 = x0(5:7);
    w0 = x0(8:10);
    v0 = x0(11:13);
    
    eR = (Rd' * R0 - R0' * Rd) / 2;
    eR = [eR(3, 2); eR(1, 3); eR(2, 1)];
    
    e0 = [p0 - pd; R0 * v0 - Rd * vd; eR; w0 - R0' * Rd * wd];
    

    [A, bmin, bmax] = gMPCConstraints(xid, e0, R0, param);
    param.P = idare(-A(1:12, 1:12), -A(1:12, param.Nx * (param.Nt + 1)+1:param.Nx * (param.Nt + 1)+param.Nu),param.Q,param.R,[], []);
    [M, q] = gMPCCost(param.Q, param.R, param.P, xid, param);
    
    solver = osqp;
    solver.setup(M, q, A, bmin, bmax,'verbose',0);
    sol = solver.solve;
    u = sol.x((param.Nt+1)*param.Nx+1:(param.Nt+1)*param.Nx+param.Nu);
    ud = coadjoint(xid(1,:)) * param.I * xid(1,:)';
    ud = ud([4,5,6,1,2,3]);
    u = u([4,5,6,1,2,3]) + ud;
end