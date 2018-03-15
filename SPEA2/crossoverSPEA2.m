function children = crossoverSPEA2(firstParent, secondParent, pc, crossoverFunction, n)
    r= rand;
    if r > pc
        children =[firstParent; secondParent];
    else
        children = crossoverFunction(firstParent,secondParent, n);
    end
end
    