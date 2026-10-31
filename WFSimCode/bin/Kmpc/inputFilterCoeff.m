function tauU = inputFilterCoeff
% inputFilterCoeff filter coefficient tau of the first-order filter of the
% turbine model, x(k+1) = (1-tau)*x(k) + tau*u(k): time constant mpc.tau = 5 s
% (MPCinitv), time step 1 s; the MPC applies the controls with this filter
% (qLPVMPCcontroller_dU_quadprog: tau = mpc.A(3,3) = 1 - tauU)
tauU = 1 - exp(-1/5);
end
