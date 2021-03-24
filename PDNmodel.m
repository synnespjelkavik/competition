% =========================================================================
% =========================================================================
% PNmodel is a function solving an ODE with
%
% res = PDNmodel(param,tRange)
%
% Input:    param :   see baseparameters()
%           tRange :  time range [days]
%        
% Output:   res :     structure containing solution of ODE
%           res.z :   spatial grid
%           res.N :   nutrient concentration (vector)
%           res.P :   phytoplankton concentration (vector)
%           res.t :   time vector
%           res.p :   parameters
% =========================================================================
% =========================================================================
function res = PDNmodel(param,tRange)
% =========================================================================
% Initial conditions
% =========================================================================
P0 = linspace(0,0,param.zb/param.dz);
P0(1:length(P0)) = 1;
D0 = linspace(0,0,param.zb/param.dz);
D0(1:length(P0)) = 1.5;
N0 = linspace(0,0,param.zb/param.dz);
N0(1:length(P0)) = 10;
y0 = [N0 D0 P0];
% =========================================================================
% Solving ODE
% =========================================================================
option = odeset('nonnegative',1:length(y0));
[t, y] = ode45(@PDNmodelDeriv, tRange, y0,option, param);
% =========================================================================
% Saving solution in structure
% =========================================================================
res.z = param.z;
res.N = y(:,1:param.nGrid);
res.D = y(:,param.nGrid+1:2*param.nGrid);
res.P = y(:,2*param.nGrid+1:3*param.nGrid);
res.t = t;
res.p = param;
% ff
% =========================================================================
end