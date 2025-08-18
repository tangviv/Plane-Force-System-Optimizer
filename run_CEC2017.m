clc
clear
close all
%%
nPop=50; % Number of population

Max_iter=500; % Maximum number of iterations


dim = 10; % Selectable 2, 10, 30, 50, 100

%%  Selection function

Function_name=3; % Function£º 1 - 30
[lb,ub,dim,fobj] = Get_Functions_cec2017(Function_name,dim);

%% Invoke the algorithm
tic
[Best_score,Best_pos,cg_curve]=PFO(nPop,Max_iter,lb,ub,dim,fobj);
toc

%% plot
figure('Position',[400 200 300 250])
semilogy(cg_curve,'Color','r','Linewidth',1)
%     plot(cg_curve,'Color','r','Linewidth',1)
title(['Convergence curve, Dim=' num2str(dim)])
xlabel('Iteration');
ylabel(['Best score F' num2str(Function_name) ]);
axis tight
grid on
box on
set(gca,'color','none')
legend('PFO')


display(['The best solution obtained by PFO is : ', num2str(Best_pos)]);
display(['The best optimal value of the objective funciton found by PFO is : ', num2str(Best_score)]);
