function [Gbest, GbestValue, ConvergenceData] = HFCLPSOLS(MaximumFEs, SwarmSize, InitPos, ModelInfor)

%-------------------- Initialisation--------------------------
me=MaximumFEs/SwarmSize;
%----------------------initialise parameters of CLPSO------------------
t = 0 : 1 / (SwarmSize - 1) : 1;
t = 5 .* t;
Pc = 0.0 + (0.5 - 0.0) .* (exp(t) - exp(t(1))) ./ (exp(t(SwarmSize)) - exp(t(1)));
m = 0 .* ones(SwarmSize, 1);
iwt = 0.9 - (1 : me) * (0.5 / me);
cc = [1.49445 1.49445]; %acceleration constants
Bound=ModelInfor.Bound;
UpperBound=Bound(:,2)';
LowerBound=Bound(:,1)';

MaximumVelocity=0.2*(UpperBound - LowerBound);
MinimumVelocity=-MaximumVelocity;
Dimension=ModelInfor.Num_WayPoints;

%w=0.9-(1:MaximumFEs).*(0.5./MaximumFEs);
ParticleSwarm=InitPos'; % SwarmSize, Dimension
Pbest=ParticleSwarm;
Velocity= MinimumVelocity + 2.*MaximumVelocity.* rand(SwarmSize, Dimension); %initialize the velocity of the particles
FitnessValue=CostFunction(ParticleSwarm', ModelInfor);
PbestValue=FitnessValue; % (1, SwarmSize)
[BestValue, GbestIndex]=min(PbestValue);
Gbest=Pbest(GbestIndex,:);  % (1, Dimension)
gbestrep = repmat(Gbest, SwarmSize, 1);
GbestValue=BestValue;
ConvergenceData=ones(1,MaximumFEs+1)*10^5;
ConvergenceData(1)=GbestValue;
fitcount=1;
stay_num = zeros(SwarmSize, 1);
ai = zeros(SwarmSize, Dimension);
f_pbest = 1 : SwarmSize;
f_pbest = repmat(f_pbest', 1, Dimension); % (SwarmSize, Dimension)

for k = 1 : SwarmSize
    ar = randperm(Dimension);
    ai(k, ar(1 : m(k))) = 1;
    fi1 = ceil(SwarmSize * rand(1, Dimension));
    fi2 = ceil(SwarmSize * rand(1, Dimension));
    fi = (PbestValue(fi1) < PbestValue(fi2)) .* fi1 + (PbestValue(fi1) >= PbestValue(fi2)) .* fi2;
    bi = ceil(rand(1, Dimension) - 1 + Pc(k));
    
    if bi == zeros(1, Dimension)
        rc = randperm(Dimension);
        bi(rc(1)) = 1;
    end
    f_pbest(k, :) = bi .* fi + (1 - bi) .* f_pbest(k, :);
end
%stop_num = 0;

%-------------------initialise parameters of LS------------------
ReductionRatio= 0.95;%reduction ratio of Quasi Entropy  减速比
PercentageParticle=1/2;% Percentage of Particles to be selected   被选择的粒子的比例
[~,sortIndex]=sort(PbestValue);    %如果 A 是向量，则 sort(A) 对向量元素进行排序   sortIndex则为PbestValue中排序后的元素在原来的位置指标
SelectedIndex=sortIndex(1:ceil(PercentageParticle*SwarmSize));  % 筛选出排前1/2的元素的位置指标
P=PbestValue(SelectedIndex)/sum(PbestValue(SelectedIndex));  % (1, SwarmSize)
InitialQuasiEntropy=-sum(P.*log(P));



for i=2: me
    
    %------------------------------------------CLPSO---------------------------------------------
    valid = [];
    for k = 1 : SwarmSize
        if stay_num(k) >= 5
            
            stay_num(k) = 0;
            ai(k, :) = zeros(1, Dimension);
            f_pbest(k, :) = k .* ones(1, Dimension);
            ar = randperm(Dimension);
            %ai(k, ar(1 : m(k))) = 1;
            fi1 = ceil(SwarmSize * rand(1, Dimension));
            fi2 = ceil(SwarmSize * rand(1, Dimension));
            fi = (PbestValue(fi1) < PbestValue(fi2)) .* fi1 + (PbestValue(fi1) >= PbestValue(fi2)) .* fi2; % find the better PbestValue between two
            bi = ceil(rand(1, Dimension) - 1 + Pc(k));
            
            if bi == zeros(1, Dimension)
                rc = randperm(Dimension);
                bi(rc(1)) = 1;
            end
            
            f_pbest(k, :) = bi .* fi + (1 - bi) .* f_pbest(k, :);
        end
        
        for dimcnt = 1 : Dimension
            pbest_f(k, dimcnt) = Pbest(f_pbest(k, dimcnt), dimcnt);
        end
        
        temp = ParticleSwarm(k, :);
        temp_ = temp + Velocity(k, :);
        
        aa(k, :) = cc(1) .* (1 - ai(k, :)) .* rand(1, Dimension) .* (pbest_f(k, :) - ParticleSwarm(k, :)) + cc(2) .* ai(k, :) .* rand(1, Dimension) .* (gbestrep(k, :) - ParticleSwarm(k, :));
        Velocity(k, :) = iwt(i) .* Velocity(k, :) + aa(k, :);
        Velocity(k, :) = (Velocity(k, :) > MaximumVelocity) .* MaximumVelocity + (Velocity(k, :) <= MaximumVelocity) .* Velocity(k, :);
        Velocity(k, :) = (Velocity(k, :) < (-MaximumVelocity)) .* (-MaximumVelocity) + (Velocity(k, :) >= (-MaximumVelocity)) .* Velocity(k, :);
        
        %% hfo
        %prev_pos=position(i,:);
        %dmax = (UpperBound - LowerBound)*sqrt(20);     % (1,20)
        rij=norm(ParticleSwarm(k, :) - Gbest)/500;
        beta0=2;
        m=2;
        gamma=1;
        beta=beta0*exp(-gamma*rij.^m);
        prev_pos1 = ParticleSwarm(k,:);
        %%
        
        ParticleSwarm(k, :) = ParticleSwarm(k, :) + Velocity(k, :);
        
      
        %%
        
        %ParticleSwarm(k, :)=(ParticleSwarm(k, :) > UpperBound).*UpperBound+ (ParticleSwarm(k, :) <= UpperBound).*ParticleSwarm(k, :);
        %ParticleSwarm(k, :)=(ParticleSwarm(k, :) < LowerBound).*LowerBound+ (ParticleSwarm(k, :) >= LowerBound).*ParticleSwarm(k, :);
        
        if (sum(ParticleSwarm(k, :) > UpperBound) + sum(ParticleSwarm(k, :) < LowerBound)) == 0
            valid = [valid k];
            fitcount = fitcount + 1;
        end
    end
    
      %% 创新点
        prev_pos2= ParticleSwarm;
       %%
    
    
    if ~isempty(valid)
        FitnessValue(valid, 1)=CostFunction(ParticleSwarm(valid, :)', ModelInfor);
        %e(valid, 1) = benchmark_func(ParticleSwarm(valid, :), problem, o, A, M, a, alpha, b);
        ConvergenceData(fitcount-length(valid)+1:fitcount)=FitnessValue(valid, 1);
        for k = 1 : length(valid)
            
            tmp = (PbestValue(valid(k)) <= FitnessValue(valid(k)));
            if tmp == 1
                stay_num(valid(k)) = stay_num(valid(k)) + 1;
            end
            temp = repmat(tmp, 1, Dimension);
            Pbest(valid(k), :) = temp .* Pbest(valid(k), :) + (1 - temp) .* ParticleSwarm(valid(k), :);
            PbestValue(valid(k)) = tmp .* PbestValue(valid(k)) + (1 - tmp) .* FitnessValue(valid(k)); %update the pbest
            if PbestValue(valid(k)) < GbestValue
                Gbest = Pbest(valid(k), :);
                GbestValue = PbestValue(valid(k));
                gbestrep = repmat(Gbest, SwarmSize, 1); %update the gbest
            end 
            
            %% 更新完Gbest和Pbest
            %hf
            ParticleSwarm(k,:)= beta.* prev_pos2(k,:) - prev_pos1;
            
            
        end
    end
    
    
   
    %%  
    
    
    %----------------------- evaluate the adaptive Local Search starting criterion --------------------------------
    [~,sortIndex]=sort(PbestValue);
    SelectedIndex=sortIndex(1:ceil(PercentageParticle*SwarmSize));
    P=PbestValue(SelectedIndex)/sum(PbestValue(SelectedIndex));  % (1, SwarmSize)
    QuasiEntropy=-sum(P.*log(P));
    
    %% 为了保证fitcount达到MaximumFEs
    
    if fitcount >= MaximumFEs  %种群乘以迭代次数
        break;
    end
    if (i == me) && (fitcount < MaximumFEs)  
        i = i-1;
    end
    %% -------------------------------------------

    %---------------------- Local search with Broyden-Fletch-Goldfarb-Shanno(BFGS)------------------------------------
    %LSMaximumFEs=min(MaximumFEs-fitcount,3000);
    LSMaximumFEs=MaximumFEs-fitcount;
    if QuasiEntropy <=ReductionRatio*InitialQuasiEntropy
        
        %[inform, x] = BFGS(Gbest',Problem,ProblemIndex,LSMaximumFEs);
        %[Gbest, GbestValue]=GreedySelection(x.p,x.f,Gbest,GbestValue);
        %inform.status
        [xopt,fopt,k, Convergence] = QuasiNewton_BFGS1(Gbest',@CostFunction,ModelInfor, LSMaximumFEs);  % set Gbest to the initial point of LS
        % 输出的是局部搜索的结果
        %k
        %options=optimset('MaxFunEvals', LSMaximumFEs);
        %[xopt,fopt] = fminsearch(@(x)Problem(x,ProblemIndex),Gbest', options);
        %fopt
        %Convergence
        %%  将局部搜索后的解与原来最好的解进行比较
        [Gbest, GbestValue]=GreedySelection(xopt,fopt,Gbest,GbestValue);
        ConvergenceData((fitcount+1):(fitcount+k))=Convergence(1:k);
        break;
        
        
    end
    
    % if fitcount >= MaximumFEs
    %     break;
    % end
    % if (i == me) && (fitcount < MaximumFEs)
    %     i = i-1;
    % end
    
end
end










