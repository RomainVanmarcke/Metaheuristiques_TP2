function children = polynomialMutation(children, pm, lowerBounds, upperBounds, n)
    r = rand(size(children));
    u = rand(size(children));
    %delta=max(upperBounds-children, children-lowerBounds);
    delta=upperBounds-lowerBounds;
    e=zeros(size(children));
    e(u<0.5) = (2*u(u<0.5)).^(1/(n+1))-1;
    e(u>=0.5) = 1-(2*(1-u(u>=0.5))).^(1/(n+1));
    children(r<pm) = children(r<pm) + delta*e(r<pm);

end