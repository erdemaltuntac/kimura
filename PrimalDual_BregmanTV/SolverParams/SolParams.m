function [Solvers] = SolParams()

tau = 50;
q_rough = .3;
MaxOuterIT_NO = 500; % Outer loop iteration number
MaxInnerIT_NO = 550; % Inner loop iteration number
lambda = 1.2; % Line search parameter

Solvers.RegPar_roughness = q_rough;
Solvers.tau = tau;
Solvers.OuterLoopNO = MaxOuterIT_NO;
Solvers.InnerLoopNO = MaxInnerIT_NO;
Solvers.lambda = lambda;



end