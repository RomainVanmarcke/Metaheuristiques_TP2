function paretoFront = main()
    T = 250; % Maximum number of generations
    t = 0;  % current generation
    N = 100; % population size
    Nb = 100; % archive size
    L = 30;  % chromosome size = number of genes
    pc = 0.9;
    pm = 1/L;
    % alpha = 0.5;
    n=20;
    %sigma = 1;
    
    M = N; % mating pool size
    
    k = round(sqrt(N+Nb));
    
    %crossoverFunction= @localArithmeticCrossover;
    crossoverFunction= @simulatedBinaryCrossover;
    
    %mutationFunction = @normalMutation;
    mutationFunction = @polynomialMutation;
    
%     feasabilityFunction = @firstFeasability;
    feasabilityFunction = @secondFeasability;
    
%     problemFunction = @SCH;
%     problemFunction = @FON;
%     problemFunction = @POL;
%     problemFunction = @KUR;
%      problemFunction = @ZDT1;
    problemFunction = @ZDT2;
%     problemFunction = @ZDT3;
%     problemFunction = @ZDT4;
%       problemFunction = @ZDT6;
        
    lowerBounds = -0; % lower bounds of the variables
    upperBounds = 1; % upper bounds of the variables
    
    pop = unifrnd(lowerBounds, upperBounds, N, L);
    archive = [];
            
   %integrer minimisation
                
    %gerer limites 
                
    while t<T
        all = [pop
            archive];
        problemValue = problemFunction(all);
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
        %disp(nNonDominated)    
        %nonDominatedDistances = distanceMatrix(:,Rmatrix==0);               
        
        
        if nNonDominated<=Nb
            [populationFitness, sortedIndices] = sort(fitnessMatrix);
            archive = all(sortedIndices(1:Nb),:);
            archiveFitness = populationFitness(1:Nb);
            paretoFront = all(Rmatrix==0,:);
                       
        else
%             distanceMatrix = distanceMatrix(Rmatrix==0,Rmatrix==0);            
            distanceMatrix = distanceMatrix(:,Rmatrix==0);            
            archive = all(Rmatrix==0,:);
            archiveFitness = fitnessMatrix(Rmatrix==0);
            %archiveValues = problemValue(Rmatrix==0,:);
            
            %j=1;
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
        
                
        %calcul du delta
%         nonDominatedValues = problemValue(Rmatrix==0,:);
%         [~,indices] = sort(nonDominatedValues(:,1));
%         if size(distanceMatrix,2)<size(nonDominatedDistances,2)
%             nonDominatedMatrix=distanceMatrix;
%             [~,nonDominatedIndices] = sort(archiveValues(:,1));
%         %c'est quoi les solutions limites?
%         %genre (1,0) et (0,1) ou celles qu'on a trouvé nous.
%         df = nonDominatedDistances(indices(1),nonDominatedIndices(1));
%         dl = nonDominatedDistances(indices(end),nonDominatedIndices(end);
%         db = 0;
%           nd = size(nonDominatedIndices,1)-1
     %      sortedDistances = zeros(nd,1);
 %         for i=1:nd
 %            sortedDistances(i)=nonDominatedDistances(nonDominatedIndices(i),nonDominatedIndices(i+1));
%             db = sortedDistances(i);
%         end
%         db = db/nd
%         delta = (df+dl+sum(abs(sortedDistances-db)))/(df+dl+nd*db)

        
        
        %Crossover
        children = zeros(M,L);
        
        for i = 1:M/2
            
            p1 = binaryTournamentSelection(archive,archiveFitness);
            p2 = binaryTournamentSelection(archive,archiveFitness);
            
            kids = crossover(p1,p2,pc, crossoverFunction, n);
            children(2*i-1,:) = kids(1,:); 
            children(2*i,:) = kids(2,:);
             
        end
%         disp(pop)
%         disp(children);
        
        %Mutation
        pop = mutationFunction(children,pm, lowerBounds, upperBounds, n);
        
        pop = feasabilityFunction(pop, lowerBounds, upperBounds);
        
%         if(t==10)
%             break
%         end
        t=t+1;   
        
        figure(1)
        plotParetoFront(paretoFront, problemFunction);   
        pause(0.01)

    end
    
            
    
end

        
       
            
            
            
            
            
            