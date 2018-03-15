function main()
    %% CONFIGURATION PART 
    % Comment/Uncomment to choose your parameters

    %PROBLEM
    problem = [];
    problem.fitnessFunction = @schaffer;
    problem = setProblemParameters(problem);
    
    %GENERAL SETTINGS
    N = 100; %Population size 
    L = problem.L; %Chromosome size
    objNumber = 2; %Number of objective functions
    Gmax = 250; %Generation max
    pc = 0.9; %Crossover probability
    pm = 1/L; %Mutation probability
    M = N; %MatingPool size
    binary = 0; %Encoding mode

    %SELECTION
    selectionFunction = @tournamentSelection; %need k
    k = 2; %size of tournament
    
    %CROSSOVER
    crossoverFunction = @simulatedBinaryCrossover;
    alpha = 0.5; %control the scope of the expansion
    
    %MUTATION
    mutationFunction = @polynomialMutation; %need n
    b = 1; %control the speed of the annealing
    sigma = 1; %standard deviation vector
    n = 20; %control parameter
    
    %FEASIBILITY
    feasibilityFunction = @firstMethod;
    
    %% EXECUTION PART
    scores = zeros(Gmax, 2*N, objNumber); %scores is a matrix of scores
    % INITIALIZATION
    pool = initialization(N, L, binary, problem.lower, problem.upper); 
    scores1 = evaluation(problem, pool, binary, objNumber);
    [~, ranks] = fastNonDominatedSort(scores1);
    matingPool= selection(selectionFunction, ranks, M, L, pool, k); %matingPool is a vector of chromosomes
    children = crossover(crossoverFunction, matingPool, pc, L, alpha); %children is a vector of chromosomes
    children = mutation(mutationFunction, children, pm, problem.lower, problem.upper, b, 1, Gmax, n, sigma);
    children = testFeasibility(feasibilityFunction, children, problem.lower, problem.upper, binary);
    
    pop = zeros(Gmax, 2*N, L);
    metrics = [];
    metrics.first = zeros(Gmax,1);
    metrics.second = zeros(Gmax,1);
    % Process all generations
    for g=1:Gmax 
        pop(g,:,:) = [pool; children];
        popg = reshape(pop(g,:, :), [2*N, L]); 
        
        % EVALUATION
        scoresg = evaluation(problem, popg, binary, objNumber);
        scores(g,:,:) = scoresg(:,:);
        [fronts, ranks] = fastNonDominatedSort(scoresg);
        
        poolIdx = [];
        f = 1;
        while (length(poolIdx) + length(fronts(f).array)) <= N
            poolIdx = [poolIdx; fronts(f).array];
            f = f + 1;
        end
        
        if (length(poolIdx) ~= N)
            scoresFrontf = [];
            for i=1:length(fronts(f).array)
                scoresFrontf = [scoresFrontf; scoresg(fronts(f).array(i),:)];
            end
            distances = crowdingDistanceAssignement(scoresFrontf, g);
            sortedLastFront = sortPartialOrder(fronts(f).array, distances);
            poolIdx = [poolIdx; sortedLastFront(1:(N - length(poolIdx)))];
        end
        
        pool = zeros(N,L);
        newRanks = zeros(N,1);
        for i=1:N
            pool(i,:) = popg(poolIdx(i),:);
            newRanks(i) = ranks(poolIdx(i));
        end
        % Variation Operators

        matingPool= selection(selectionFunction, newRanks, M, L, pool, k); %matingPool is a vector of chromosomes
        children = crossover(crossoverFunction, matingPool, pc, L, alpha); %children is a vector of chromosomes
        children = mutation(mutationFunction, children, pm, problem.lower, problem.upper, b, g, Gmax, n, sigma);
        children = testFeasibility(feasibilityFunction, children, problem.lower, problem.upper, binary);

        % Find pareto front and distance metric foreach gen
        firstFrontIdx = fronts(1).array;
        scoresFirstFront = zeros(length(firstFrontIdx),2);
        for i=1:length(firstFrontIdx)
            scoresFirstFront(i,:) = scoresg(firstFrontIdx(i),:);
        end
        
        scoresPool = evaluation(problem, pool, binary, objNumber);
        if length(problem.pareto) > 1
            scoresPareto= evaluation(problem, problem.pareto, binary, objNumber);
            metrics = computeMetrics(metrics, scoresPareto, scoresPool, scoresFirstFront, g);
        end
        
        
    end
    
    %% DISPLAY PART
    fprintf('Result for Gen %d, Popsize %d \n', g,N);
%     scoresLastGen = evaluation(problem, pool, binary, objNumber);
    displayResult(problem, scoresFirstFront, scoresPareto, metrics);
    
end
    

