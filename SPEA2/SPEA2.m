function [paretoFrontAllG, scores ] = SPEA2(problem)
    T = problem.Gmax; % Maximum number of generations
    t = 0;  % current generation
    N = problem.N; % population size
    Nb = 100; % archive size
    L = problem.L;  % chromosome size = number of genes
    pc = 0.9;
    pm = 1/L;
    % alpha = 0.5;
    n=20;
    %sigma = 1;
    
    M = N; % mating pool size
    
    k = round(sqrt(N+Nb));

    crossoverFunction= @simulatedBinaryCrossover;
    
    mutationFunction = @polynomialMutationSPEA2;
    
%     feasabilityFunction = @firstFeasability;
    feasabilityFunction = @secondFeasability;
    
    problemFunction = problem.alias;
        
    lowerBounds = problem.lower; % lower bounds of the variables
    upperBounds = problem.upper; % upper bounds of the variables
    
    if length(lowerBounds) == 1
        lowerBounds = [lowerBounds, lowerBounds];
        upperBounds = [upperBounds, upperBounds];
    end
        
    pop = unifrnd(lowerBounds(2),upperBounds(2),N,L);
    pop(:,1) = unifrnd(lowerBounds(1),upperBounds(1),N,1);

    archive = [];

    scores(1).array = []; 
    paretoFrontAllG(1).array = [];
    while t<T+1
        all = [pop
            archive];
        problemValue = problemFunction(all);
        if t > 0
            scores(t).array = problemValue(:,:);
        end
        nall = size(all,1);
        dominationMatrix = false(nall, nall);
        
        Smatrix = zeros(nall, 1);
        
        for i=1:nall
            for j=i+1:nall
                
                if dominates(problemValue(i,:),problemValue(j,:))
                    Smatrix(i)= Smatrix(i)+1;
                    dominationMatrix(i,j) = true;
                    
                elseif dominates(problemValue(j,:),problemValue(i,:))
                    Smatrix(j) = Smatrix(j)+1;
                    dominationMatrix(j,i)=true;
                    
                end
                
            end
        end
        
            
        Rmatrix = zeros(nall, 1);
        for i = 1:nall
            Rmatrix(i) = sum(Smatrix(dominationMatrix(:,i)));

        end
        
        distanceMatrix = pdist2(problemValue, problemValue, 'seuclidean');
        distanceMatrix = sort(distanceMatrix);

        sigmaKmatrix = distanceMatrix(k,:);
        
        Dmatrix = 1./(sigmaKmatrix'+2);
        fitnessMatrix = Rmatrix + Dmatrix;
        
        nNonDominated = sum(Rmatrix==0);    
        
        if nNonDominated<=Nb
            [populationFitness, sortedIndices] = sort(fitnessMatrix);
            archive = all(sortedIndices(1:Nb),:);
            archiveFitness = populationFitness(1:Nb);
            paretoFront = all(Rmatrix==0,:);
                       
        else   
            distanceMatrix = distanceMatrix(:,Rmatrix==0);            
            archive = all(Rmatrix==0,:);
            archiveFitness = fitnessMatrix(Rmatrix==0);

            j=2;
            while size(archive,1)>Nb
                while min(distanceMatrix(j,:)) == max(distanceMatrix(j,:)) && j<size(distanceMatrix,1)
                    j=j+1;
                end
                
                [~,i] = min(distanceMatrix(j,:));
                
                archive(i,:)=[];
                archiveFitness(i)=[];
                %archiveValues(i,:)=[];
                distanceMatrix(:,i) = [];
                
            end
            paretoFront = archive;
            
        end   
        
        %Crossover
        children = zeros(M,L);
        
        for i = 1:M/2
            
            p1 = binaryTournamentSelection(archive,archiveFitness);
            p2 = binaryTournamentSelection(archive,archiveFitness);
            
            kids = crossoverSPEA2(p1,p2,pc, crossoverFunction, n);
            children(2*i-1,:) = kids(1,:); 
            children(2*i,:) = kids(2,:);
             
        end
        
        %Mutation
        pop = mutationFunction(children,pm, lowerBounds, upperBounds, n);
        
        pop = feasabilityFunction(pop, lowerBounds, upperBounds);

        if t > 0 
            paretoFrontAllG(t).array = problemFunction(paretoFront);
        end 
        t=t+1;   
        
%         figure(1)
%         plotParetoFront(paretoFront, problemFunction);   
%         pause(0.01)

    end
    
            
    
end

        
       
            
            
            
            
            
            
