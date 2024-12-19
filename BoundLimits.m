function X=BoundLimits(X,lower,upper)

    X(X>upper)=upper;
    X(X<lower)=lower;
end