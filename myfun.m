function [f,g] = myfun(x,CostFunction, ModelInfor) 
%  x : initial point  (D,1)
%  f : function value
%  g : gradient
%
f=CostFunction(x, ModelInfor);    
%f= feval(Problem, x , ProblemIndex);
g = gradFunction(CostFunction, ModelInfor, x);
    
end
