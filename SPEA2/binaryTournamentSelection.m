function selectedIndividual = binaryTournamentSelection(archive, archiveFitness)
    
    r = randsample(size(archive,1),2);
    
    firstIndividual = r(1);
    secondIndividual = r(2);
    
    if archiveFitness(firstIndividual)<archiveFitness(secondIndividual)
        selectedIndividual = archive(firstIndividual,:);
    else
        selectedIndividual = archive(secondIndividual,:);
    end
    
end
    