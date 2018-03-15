function pop = secondFeasability(pop, lowerBounds, upperBounds)
        pop(pop<lowerBounds) = 2*lowerBounds-pop(pop<lowerBounds);
        pop(pop>upperBounds) = 2*upperBounds-pop(pop>upperBounds);
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
