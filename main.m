clc;
clear;
addpath("../improved_CTCR/");

need_reference = input('need_reference = ');

P = input('P = ');
D = input('D = ');

graph_start = input('start_points = ');
graph_end = input('end_points = ');
Graph = graph(graph_start, graph_end);

% build graph and the whole system
% need_reference = 0;
% P = 5;
% D = 1;
% start = [1,1,2,2,3,3,4,2];
% ends = [2,3,4,5,4,5,5,3];
% Graph = graph(start, ends);


n = size(Graph.Nodes);
n = n(1);
Delta = diag(Graph.degree);

A_Gamma = adjacency(Graph);
C = Delta\A_Gamma;
[T, Lambda] = eig(C);

D_max = zeros(1,n);
for i = 1:n
    if Lambda(i,i) <= 0.99999 && Lambda(i,i) >= -0.99999
        D_max(i) = sqrt(2*P*(1-sqrt(1-Lambda(i,i)^2))/(1-Lambda(i,i)^2));
    elseif Lambda(i,i) == 1
        D_max(i) = sqrt(2*P*(1-sqrt(1-(Lambda(i,i)-0.0001)^2))/(1-(Lambda(i,i)-0.0001)^2));
    else
        D_max(i) = sqrt(2*P*(1-sqrt(1-(Lambda(i,i)+0.0001)^2))/(1-(Lambda(i,i)+0.0001)^2));
    end
end

if need_reference == 1
    plot_tau_D(Lambda, D_max,n,P);
end

exi_lambda = min(eig(Lambda));
% count the system's stable tau zone 
if exi_lambda <= -0.9999
    [exi_result_table, ~] = CTCR([0,1;-P,-D], [0,0;P*(exi_lambda+0.0001), D*(exi_lambda+0.0001)]);
else
    [exi_result_table, ~] = CTCR([0,1;-P,-D], [0,0;P*exi_lambda, D*exi_lambda]);
end
[max_result_table, ~] = CTCR([0,1;-P,-D], [0,0;P*(1-0.0001),D*(1-0.0001)]);
stable_max_tau = min(exi_result_table(1,1), max_result_table(1,1));

% fast reaching consensus tau
fastest_tau = 0;
rightmost_root = 0;
for tau = 0: stable_max_tau/10: stable_max_tau
    rightmost_roots_for_lambda = zeros(1,n);
    for k = 1:n
        lambda = Lambda(k,k);
        if lambda == 1
            lambda = 1-0.0001;
        end
        if tau == 0
            rightmost_roots_for_lambda(k) = max(real(QPmR([-10,0,-100,100],[1,D,P;0,-D*lambda,-P*lambda],[0;tau+0.001])));
        else
            rightmost_roots_for_lambda(k) = max(real(QPmR([-10,0,-100,100],[1,D,P;0,-D*lambda,-P*lambda],[0;tau])));
        end
    end
    if tau == 0
        rightmost_root = max(rightmost_roots_for_lambda);
    else
        if max(rightmost_roots_for_lambda) < rightmost_root
            rightmost_root = max(rightmost_roots_for_lambda);
            fastest_tau = tau;
        end
    end
end

for tau = fastest_tau - stable_max_tau/10: stable_max_tau/50: fastest_tau + stable_max_tau/10
    rightmost_roots_for_lambda = zeros(1,n);
    for k = 1:n
        lambda = Lambda(k,k);
        if lambda == 1
            lambda = 1-0.0001;
        end
        if tau == 0
            rightmost_roots_for_lambda(k) = max(real(QPmR([-10,0,-100,100],[1,D,P;0,-D*lambda,-P*lambda],[0;tau+0.001])));
        else
            rightmost_roots_for_lambda(k) = max(real(QPmR([-10,0,-100,100],[1,D,P;0,-D*lambda,-P*lambda],[0;tau])));
        end
    end
    if tau == 0
        rightmost_root = max(rightmost_roots_for_lambda);
    else
        if max(rightmost_roots_for_lambda) < rightmost_root
            rightmost_root = max(rightmost_roots_for_lambda);
            fastest_tau = tau;
        end
    end
end

for tau = fastest_tau - stable_max_tau/50: stable_max_tau/250: fastest_tau + stable_max_tau/50
    rightmost_roots_for_lambda = zeros(1,n);
    for k = 1:n
        lambda = Lambda(k,k);
        if lambda == 1
            lambda = 1-0.0001;
        end
        if tau == 0
            rightmost_roots_for_lambda(k) = max(real(QPmR([-10,0,-100,100],[1,D,P;0,-D*lambda,-P*lambda],[0;tau+0.001])));
        else
            rightmost_roots_for_lambda(k) = max(real(QPmR([-10,0,-100,100],[1,D,P;0,-D*lambda,-P*lambda],[0;tau])));
        end
    end
    if tau == 0
        rightmost_root = max(rightmost_roots_for_lambda);
    else
        if max(rightmost_roots_for_lambda) < rightmost_root
            rightmost_root = max(rightmost_roots_for_lambda);
            fastest_tau = tau;
        end
    end
end

disp('for given P,D, the stable tau zone is :');
disp(['0-', num2str(stable_max_tau)]);
disp('to have a fastest reaching consensus speed, you can set tau to be:')
disp(num2str(fastest_tau));
disp('rightmost root:');
disp(num2str(rightmost_root));
