function Ur = urFromPower(sol, Wp)
% urFromPower rotor effective wind speeds estimated from the measured
% power, the thrust coefficient CT' and the yaw angle (WFSim power model,
% Actuator.m): P = powerscale(phi) * 0.5 * rho * A * CT' * mean(Ue^3), so
% Ur = (P / (powerscale(phi) * 0.5 * rho * A * CT'))^(1/3)
%
% As the power uses the mean of Ue^3 over the rotor, the estimate is the
% cubic mean and thus larger than the measured Ur = mean(Ue) for a
% non-uniform inflow (partial wake of the downstream turbine)
%
% Inputs
% - sol: WFSim solution with sol.turbine.power, .CT_prime, .Phi
% - Wp: WFSim parameters (Wp.site.Rho, Wp.turbine.Drotor, .powerscale)
%
% Output
% - Ur: estimated rotor effective wind speeds (N x 1)

cpVec = [0,15,30,45; 0.5, 0.48,0.4,0.26]; % as Actuator.m
cpYaw = interp1(cpVec(1,:), cpVec(2,:)/cpVec(2,1), abs(sol.turbine.Phi(:)));
Ar = pi*(0.5*Wp.turbine.Drotor)^2;
cp = Wp.turbine.powerscale(:).*cpYaw*0.5*Wp.site.Rho*Ar;
Ur = nthroot(max(sol.turbine.power(:),0)./(cp.*sol.turbine.CT_prime(:)), 3);
