clc
clear
close all
%%
nPop=50; % Number of population
Max_iter=500; % Maximum number of iterations
dim = 10; % Selectable 2, 10, 20

%%  Selection function
Function_name=11; % Function£º 1 - 12
[lb,ub,dim,fobj] = Get_Functions_cec2022(Function_name,dim);

%% Invoke the algorithm
Optimal_results={}; % Optimal results
index = 1;
% PFO
tic
[Best_score,Best_pos,cg_curve]=PFO(nPop,Max_iter,lb,ub,dim,fobj);
Optimal_results{1,index}="PFO";
Optimal_results{2,index}=cg_curve;
Optimal_results{3,index}=Best_score;
Optimal_results{4,index}=Best_pos;
Optimal_results{5,index}=toc;
index = index +1;


%% plot
figure
semilogy(cg_curve,'Color','r','Linewidth',1)
% for i = 1:size(Optimal_results, 2)
%     plot(Optimal_results{2, i},'Linewidth',2)
%     hold on
% end
title(['Convergence curve, Dim=' num2str(dim)])
xlabel('Iteration');
ylabel(['Best score F' num2str(Function_name) ]);
axis tight
grid on
box on
set(gcf,'Position',[400 200 400 250])
legend(Optimal_results{1, :})

