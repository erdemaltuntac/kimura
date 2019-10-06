 function [u_rec, MfitErr_val,RelatErr_val] = ... 
     PrimalDual_BregmanTV(T, u_0 , point_tr , v_delta , delta , sigma , Geometry)
%%                      PrimalDual_BregmanTV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Erdem Altuntac
%
%                       Universite Libre de Bruxelles
%                       Department of Mathematics
%
%                       e-mail: erdem.altuntac@ulb.ac.be
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bregman iterated nested primal-dual algorithm;
% -------------------------------------------------------------------------
% Bregman distance associated with non-smooth TV.
% Final reconstruction is obtained by a simple line-search.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SOLVER PARAMETERS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Solvers] = SolParams();
tau = Solvers.tau;
MaxIT_NO = Solvers.OuterLoopNO;
InnerLoop_NO = Solvers.InnerLoopNO;
lambda = Solvers.lambda;
mu = .5/(sigma^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ASSIGNMENTS AND VECTORIZATION
% State/Pre-image Space
u_init = u_0;
point_tr = point_tr(:);
T_T = T';
%--------------------------------------------------------------------------
% Activate below lines in case you need to test with different amount of 
% noise input
%--------------------------------------------------------------------------
% Measurement Space
%v_truth = problem.Measurement.y_tr;
%[v_delta , delta] = addnoise(v_truth);
%--------------------------------------------------------------------------
%% OUTER LOOP INITIATION AND FIRST ERROR ANALYSIS
Out_ItStep = 1;
DiscRad = tau * delta;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Mfit0 , Mfit_norm] = Misfit(T, u_0, v_delta);
[Grad , Grad_T , Dx , Dy] = FiniteDiff2D(Geometry);
[Grad_u0, u0_x , u0_y , TV_u0] = DiscrDiff_Val(u_0 , Dx , Dy);
[w_00, w00_x , w00_y] = SubDiff_L1(u0_x, u0_y);
[w_10 , w10_x , w10_y] = VecAssign(w00_x , w00_y);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Transpose of the gradient %%%%%%%%%%%%%%%%%%%%%%
DxT = -Dx';
DyT = -Dy';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Subdifferential of the Bregman distance %%%%%%%%%%%%
[DT_w00] = GradT_Dual(w00_x , w00_y , DxT , DyT);
[DT_w10] = GradT_Dual(w10_x , w10_y , DxT , DyT);
DT_wDiff10 = DT_w10 - DT_w00;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Initial diagnostics %%%%%%%%%%%%%%%%%%%%%%%

MfitErr_val(Out_ItStep) = Mfit_norm;
RelatErr_val(Out_ItStep) = norm(u_0 - point_tr)/norm(point_tr);

AA = zeros(MaxIT_NO,1);
BB = zeros(MaxIT_NO,1);
CC = zeros(MaxIT_NO,1);
DD = zeros(MaxIT_NO*InnerLoop_NO,1);
EE = zeros(MaxIT_NO,1);
GG = zeros(MaxIT_NO*InnerLoop_NO,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




while(Mfit_norm >= DiscRad && Out_ItStep < MaxIT_NO)
    %% BEGIN
    
    Out_ItStep = Out_ItStep + 1;
    OuterIterationNO_Info = " >>>>>>>>>>>>>>>>>>>>>>>>> Current outer loop iteration step is: %d ";
    sprintf(OuterIterationNO_Info,Out_ItStep)
    
    %% GRADIENT OF THE DATA-FIT TERM GETS UPDATED
    
    Grad_Mfit = T_T * Mfit0;
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% REGULARIZATION PARAMETER UPDATE
    
    alpha_1 = (delta)/(Out_ItStep^2);
    regpar = alpha_1;
           
    disp(['Updated regpar is ',num2str(regpar)])
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%% Inner Loop: Dual Update %%%%%%%%%%
    Inn_ItStep = 0;
    
    while(Inn_ItStep < InnerLoop_NO)
        
        Inn_ItStep = Inn_ItStep + 1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp(['Current inner loop iteration step is >> ',num2str(Inn_ItStep)])
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%                     PRIMAL VARIABLE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [u_01] = PrimalUpdate(u_0 , Grad_Mfit , regpar , mu , DT_wDiff10);
        u_primal = u_01;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%        Some error analysis in advance before the recurrance
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        AA(Out_ItStep) = RelatErr(u_primal , u_0);
        disp(['Relative iterative error; u_{i}^{j} against u_{i} -------> ',num2str(AA(Out_ItStep))])
                
        FG = RelatErr(u_primal , point_tr);
        disp(['Relative error value u_{i}^{j} against u^{+} -------> ',num2str(FG)])
        GG(Inn_ItStep*Out_ItStep) = FG;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%                       DUAL VARIABLE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        [Grad_uPrime, uPrime_x , uPrime_y , TV_uPrime] = DiscrDiff_Val(u_primal , Dx , Dy);
        nu = (delta/Out_ItStep)^Out_ItStep;
        [w_11 , w11_x , w11_y] = DualUpdate(w10_x , w10_y , nu , regpar , uPrime_x , uPrime_y);
        
        ED = RelatErr(w_11,w_10);
        DD(Out_ItStep*Inn_ItStep) = ED;
        disp(['Relative error value w_{i}^{j+1} against w_{i}^{j} -------> ',num2str(ED)])
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Recursion of w and calculating the subdiff of Bregman
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [w_10 , w10_x , w10_y] = VecAssign(w11_x , w11_y);
        [DT_w10] = GradT_Dual(w10_x , w10_y , DxT , DyT);
        DT_wDiff10 = DT_w10 - DT_w00;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%% Line Search %%%%%%%%%%%%%
    u_rec = LineSearch(u_0, u_01, lambda);
    
    %%%%%%%%%%% Display %%%%%%%%%%%%%%%%%%%%%
    
%     filename = 'Star.gif';
%     imagesc(reshape(u_rec,m_x,n_x),[0,max(x_tr(:))]), colorbar, colormap(pink)
%     title('Reconstruction')
%     drawnow
%     frame = getframe;
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,512);
%     
%     if Out_ItStep == 2
%         imwrite(imind,cm,filename,'gif','DelayTime',0,'LoopCount',inf);
%     else
%         imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0);
%     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% Error values on the pre-image space %%%%%%%%%
    
    ZZ = RelatErr(u_rec , point_tr);
    disp(['>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Relative error value; u_{i+1} against u^{+} -------> ',num2str(ZZ)])
    RelatErr_val(Out_ItStep) = ZZ;
    
    AB = SuccErr(u_rec , u_01);
    disp(['Successive error value; u_{i+1} against u_{i}^{j} -------> ',num2str(AB)])
    BB(Out_ItStep) = AB;
    
    BC = SuccErr(u_rec , u_init);
    disp(['Successive error value; u_{i+1} against u^{0} -------> ',num2str(BC)])
    CC(Out_ItStep) = BC;
    
    CE = RelatErr(u_rec , u_0);
    disp(['Relative error value; u_{i+1} against u^{i} -------> ',num2str(CE)])
    EE(Out_ItStep) = CE;
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%% Update on the data-fit %%%%%%%%%%%%%%%%
    [Mfit0 , Mfit_norm] = Misfit(T, u_rec, v_delta);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Recursion of u and w 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    u_0 = u_rec;
    [Grad_u0 , u0_x , u0_y , TV_u0] = DiscrDiff_Val(u_0, Dx , Dy);
    [w_00, w00_x , w00_y] = SubDiff_L1(u0_x , u0_y);
    
    [DT_w00] = GradT_Dual(w00_x , w00_y , DxT , DyT);
    DT_wDiff10 = DT_w10 - DT_w00;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%% Discrepancy on the image space %%%%%%%%%%%%%%
    
    MfitErr_val(Out_ItStep) = Mfit_norm;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
end


end




