function Return = gradFunction(CostFunction, ModelInfor, seed)
%
% numerical computation of gradient
% this allows automatic gradient computation
% 
% seed: the seed (Dimension,1)
% first forward finite difference
% hstep = 0.001; - programmed in
%
hstep = 0.0000001;
Dimension = length(seed);
f= CostFunction(seed, ModelInfor); 
for i = 1:Dimension
   xs = seed;
   xs(i) = xs(i) + hstep;
   %CostFunction(seed, ModelInfor)
   gradx(i)= (CostFunction(xs, ModelInfor) -f)/hstep;
end
Return = gradx'; % (Dimension,1)