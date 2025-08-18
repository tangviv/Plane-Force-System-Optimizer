
clear all 
clc

SearchAgents_no=30; % Number of search agents

Function_name='F8'; % Name of the test function that can be from F1 to F23

Max_iteration=500; % Maximum numbef of iterations

% Load details of the selected benchmark function
[lb,ub,dim,fobj]=Get_Functions_details(Function_name);

[Best_score,Best_pos,PFO_cg_curve]=PFO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);

figure('Position',[269   240   660   290])
%Draw search space
subplot(1,2,1);
func_plot(Function_name);
title('Parameter space')
xlabel('x_1');
ylabel('x_2');
zlabel([Function_name,'( x_1 , x_2 )'])

%Draw objective space
subplot(1,2,2);
semilogy(PFO_cg_curve,'Color','r')
title('Objective space')
xlabel('Iteration');
ylabel('Best score obtained so far');

axis tight
grid on
box on
legend('PFO')

display(['The best solution obtained by PFO is : ', num2str(Best_pos)]);
display(['The best optimal value of the objective funciton found by PFO is : ', num2str(Best_score)]);

        



