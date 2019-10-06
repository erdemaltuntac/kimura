
%%                      Solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Erdem Altuntac
%
%                       Universite Libre de Bruxelles
%                       Department of Mathematics
%
%                       e-mail: erdem.altuntac@ulb.ac.be
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% LOAD THE PROBLEM
clear;close all;clc
load('Star_65536X65536_noise1');

T = problem.Measurement.ForwardOp;
[U,S,V] = svds(T);
sigma = max(max(S));

clearvars -except problem sigma
close all

x_tr = problem.Model.x_tr;
[m_x,n_x] = size(x_tr);
y_noise = problem.Measurement.MeasNoisy(:);
delta = problem.Measurement.delta;
Geometry = problem.Geometry;

%% ASSIGNMENTS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T = problem.Measurement.ForwardOp;
[Row_T , Col_T] = size(T);
u_0 = zeros(Col_T,1);
u_init = u_0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SOLVERS

[u_rec, MfitErr_val,RelatErr_val] = ... 
    PrimalDual_BregmanTV(T, u_0 , x_tr , y_noise , delta , sigma , Geometry);

%% EROR ANALYSIS PROFILES

figure(4)
subplot(1,3,1)
imagesc(x_tr), colorbar
title('The original image - Ground truth')

subplot(1,3,2)
imagesc(y_noise), colormap(bone), colorbar
title('Sinogram')

subplot(1,3,3)
imagesc(reshape(u_rec,m_x,n_x),[0,max(x_tr(:))]), colorbar, colormap(gray)
title('Reconstruction from noisy sinogram measurement')

figure(5)
title('Diagnostics; Convergence & Error Analysis')

subplot(2,3,1)
semilogy(MfitErr_val(1:length(MfitErr_val)),'-p')
grid on 
hold on
title('Discrepancy; Misfit')
xlabel('Number of iterations')
ylabel('||T u_{i+1} - v||')

subplot(2,3,2)
semilogy(RelatErr_val(1:length(RelatErr_val)),'-p')
grid on
hold on
title('Relative Error Value of the Reconstruction; u_{rec} and u^{+}')
xlabel('Number of iterations')
ylabel('||u_{i+1} - u^{+}||/|| u^{+} ||')

subplot(2,3,3)
imagesc(u_init), colormap(gray), colorbar
title('Initial guess')

subplot(2,3,4)
imagesc(x_tr), colormap(gray), colorbar
title('True data on the pre-image space')

subplot(2,3,5)
imagesc(reshape(u_rec,m_x,n_x),[0,max(x_tr(:))]), colorbar, colormap(gray)
title('Reconstruction')


%% SAVE THE RESULTS IF NECESSARY

Results.Recon = u_rec;
Results.MfitErr_val = MfitErr_val;
Results.RelatErr_val = RelatErr_val;
Results.InitGuess = u_init;
Results.problem = problem;


