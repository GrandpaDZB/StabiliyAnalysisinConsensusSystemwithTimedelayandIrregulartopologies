function [] = plot_tau_D(Lambda,D_max,n,P)
%plot the relationship between tau and D
%   need eigenvalue matrix, D_max vector, number of the agents, param P
% plot

ft = fittype(@(a,b,x)a*exp(b*x));

figure

for j = 1:n
lambda = Lambda(j,j);
sample_number = 5;
plot_limit = 5;

sample_D = zeros(1,sample_number);
sample_t = zeros(2,sample_number);

if lambda >= 0.99999 || lambda <= -0.99999
    for i = 1:sample_number
        if i == 1
            sample_D(i) = plot_limit/(sample_number-1)*(i-1)+0.0001;
        else
            sample_D(i) = plot_limit/(sample_number-1)*(i-1);
        end
        A = [0,1;-P,-sample_D(i)];
        if lambda >= 0.99999
            B = [0,0;(lambda-0.00001)*P,(lambda-0.00001)*sample_D(i)];
        else
            B = [0,0;(lambda+0.00001)*P,(lambda+0.00001)*sample_D(i)];
        end    
        [result_table, ~] = CTCR(A,B);

        sample_t(2,i) = result_table(1,1);
        sample_t(1,i) = result_table(2,1);
    end
else
    for i = 1:sample_number
        if i == 1
            sample_D(i) = D_max(j)/(sample_number-1)*(i-1)+0.0001;
        else
            sample_D(i) = D_max(j)/(sample_number-1)*(i-1)-0.001;
        end
        A = [0,1;-P,-sample_D(i)];
        B = [0,0;(lambda)*P,(lambda)*sample_D(i)];
        [result_table, ~] = CTCR(A,B);

        sample_t(2,i) = result_table(1,1);
        sample_t(1,i) = result_table(2,1);
    end
end

if D_max(j) < plot_limit
    plot_x = (0:D_max(j)/50:D_max(j));
else
    plot_x = (0:plot_limit/50:plot_limit);
end
disp(sample_t);
fit_params = polyfit(sample_D, sample_t(1,:),4);
plot_y1 = fit_params(1)*plot_x.*plot_x.*plot_x.*plot_x+fit_params(2)*plot_x.*plot_x.*plot_x+fit_params(3).*plot_x.*plot_x+fit_params(4).*plot_x+fit_params(5);
fit_params = polyfit(sample_D, sample_t(2,:),4);
plot_y2 = fit_params(1)*plot_x.*plot_x.*plot_x.*plot_x+fit_params(2)*plot_x.*plot_x.*plot_x+fit_params(3).*plot_x.*plot_x+fit_params(4).*plot_x+fit_params(5);

if lambda >= 0.99999 || lambda <= -0.99999
    plot(plot_x, plot_y2,'-r','LineWidth', 2);
    text(plot_limit, plot_y2(end), ['l=',num2str(lambda)])
else
    plot(plot_x, plot_y1,'-b',plot_x, plot_y2,'-r','LineWidth', 2);
    text(D_max(j),(plot_y1(end)+plot_y2(end))/2, ['Dmax  ','l=',num2str(lambda)]);
end
xlabel('D');
ylabel('\tau /s');
hold on;
end
end

