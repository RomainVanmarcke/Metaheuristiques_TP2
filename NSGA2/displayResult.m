function displayResult(problem, scoresFirstFront, scoresPareto, metrics)

    % Display all chromosomes of the first front for the last generation
    figure;
    grid on
    subplot(2,2,[1,2])
    plot(scoresFirstFront(:,1),scoresFirstFront(:,2),'ro');

    if length(problem.pareto) > 1
        hold on;
        line(scoresPareto(:,1),scoresPareto(:,2));
        hold off;
        xlabel('f_1');
        ylabel('f_2');
        legend('NSGA-II Pareto Front', 'Optimal Pareto Front');
        title(['Pareto Front for ', func2str(problem.fitnessFunction), ' problem']);


        subplot(2,2,3)
        X = 1:length(metrics.first);
        plot(X,log10(metrics.first(X)));
        xlabel('Generation Number');
        ylabel('Average Distance (log10)');
        legend('NSGA-II Average Distance');
        title(['Distance metric for ', func2str(problem.fitnessFunction)]);
        
        subplot(2,2,4)
        X = 1:length(metrics.second);
        plot(X,(metrics.second(X)));
        xlabel('Generation Number');
        ylabel('Delta');
        legend('NSGA-II Diversity metric');
        title(['Diversity metric for ', func2str(problem.fitnessFunction)]);
    end
end