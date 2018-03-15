function pop = secondFeasability(pop, lowerBounds, upperBounds)
        pop1=pop(:,1);
        pop1(pop1<lowerBounds(1)) = 2*lowerBounds(1)-pop1(pop1<lowerBounds(1));
        pop1(pop1>upperBounds(1)) = 2*upperBounds(1)-pop1(pop1>upperBounds(1));

        popn=pop(:,2:end);
        popn(popn<lowerBounds(2)) = 2*lowerBounds(2)-popn(popn<lowerBounds(2));
        popn(popn>upperBounds(2)) = 2*upperBounds(2)-popn(popn>upperBounds(2));

        pop = [pop1 popn];
%         
%         a =zeros(size(pop));
%         
%         for i = 1:size(pop,1)
%             for j = 1:size(pop,2)
%                 if pop(i,j)<lowerBounds
%                     a(i,j) = 2*lowerBounds-pop(i,j);
%                 elseif pop(i,j)>upperBounds
%                     a(i,j) = 2*upperBounds-pop(i,j);
%                 else 
%                     a(i,j)=pop(i,j);
%                 end
%             end
%         end
% %                     
%         disp(a)
%         disp(pop)
end
