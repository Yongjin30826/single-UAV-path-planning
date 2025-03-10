% Revised: 12.2024
% This is HFCLPSOLS to optimize problems of UAV path planning, based on the paper 
% "Yongjin Wang, Shuang Geng, Ben Niu. Hybrid Algorithm Based on Comprehensive Learning Particle Swarm Optimisation with 
% Local Search and Firefly Algorithm for UAV path planning.IEEE,CSIS-IAC(2024)"
clear;
clc;
% AlgorithmName={'CLPSOLSBFGS','PSO','CLPSOLSBFGS3'};
AlgorithmName={'HFCLPSOLS'};
% algorithm parameters  算法的参数
TrialTimes=1;   % 训练次数
MaximumFEs=15000;   %迭代次数
%MaximumFEs=50000;   %迭代次数
SwarmSize=40;      % 种群数量


%Dimension=20;   % the number of waypoints   垂直于起始点与目标点直线的数量

TimeUsed=zeros(1, TrialTimes);
ConvergenceData=zeros(TrialTimes,MaximumFEs+1);  %record the convergence data for each run  记录每次收敛的数据
for Dimension=[20 ]
BestPathAndCost=zeros(Dimension+1,TrialTimes);   %record the best fitness value for each run    记录每次飞行的最好的适应度成本值    
file_result=strcat('Iter',int2str(MaximumFEs),'Size',int2str(SwarmSize),'Dim',int2str(Dimension),'Result.txt'); % record results of  all functions by all algorithms (average, std, time) 记录所有算法所有函数的结果
find_file_result=fopen(file_result,'a+');  % fopen（filename,permisison）打开或创建要写入的新文件。追加数据到文件末尾。

for TaskIndex=[1]  %1：5
    [TaskInfor, ThreatInfor, ObstacleInfor ]=EnvironmentInfor(TaskIndex);
    ModelInfor=CordinateTransformation(TaskInfor, ThreatInfor, ObstacleInfor, Dimension); % coordinate transformaion, let start point be (0,0) 坐标变换
    % Initialise positions
    Bound=ModelInfor.Bound;
    
    for AlgorithmIndex=[1]  %
        file_convergence_data=strcat(char(AlgorithmName(AlgorithmIndex)),'Prob',int2str(TaskIndex),'Dim',int2str(Dimension),'Data.txt'); %record the convergence data  记录收敛数据
        find_file_convergence_data=fopen(file_convergence_data,'w');  %打开或创建要写入的新文件，并放弃现有的内容
        file_eachTrialResult_data=strcat(char(AlgorithmName(AlgorithmIndex)),'Prob',int2str(TaskIndex),'Dim',int2str(Dimension),'Result.txt'); %record the best position and fitness values found by each trial  记录每次训练的最好的粒子和适应度值
        find_file_eachTrialResult_data=fopen( file_eachTrialResult_data,'w');
        file_bestpath=strcat(char(AlgorithmName(AlgorithmIndex)), 'Prob',int2str(TaskIndex), 'Dim',int2str(Dimension),'Path.txt'); % record the best path among all the trials by all algorithms  记录所有算法中在所有训练次数中最好的路径
        find_file_bestpath=fopen(file_bestpath,'w');
        
        Algorithm=str2func(char(AlgorithmName(AlgorithmIndex)));
        for TrialIndex=1:TrialTimes
            TimeStart=tic; % record the running time  记录运行时间
            InitPos=Bound(:,1)+rand(Dimension, SwarmSize).*(Bound(:,2)-Bound(:,1));
            [Gbest, GbestValue, GbestHistory]=feval(Algorithm, MaximumFEs, SwarmSize, InitPos, ModelInfor);
            BestPathAndCost(1, TrialIndex)=GbestValue; % Optimal cost
            
            BestPathAndCost(2:Dimension+1, TrialIndex)= Gbest; % 
            for i=1:MaximumFEs
                if GbestHistory(i)<GbestHistory(i+1)
                    GbestHistory(i+1)=GbestHistory(i);
                end   
            end
            ConvergenceData(TrialIndex,:)=GbestHistory(1:MaximumFEs+1);
            TimeUsed(TrialIndex)=toc(TimeStart);
        end
        AverageTime=mean(TimeUsed);
        StdTime=std(TimeUsed);
        MeanOfBestCost=mean(BestPathAndCost(1,:));
        StdOfBestCost=std(BestPathAndCost(1,:));
        MaxOfBestCost=max(BestPathAndCost(1,:));
        [MinOfBestCost, BestIndex]=min(BestPathAndCost(1,:));
        BestY=BestPathAndCost(2:Dimension+1, BestIndex);
        BBestPath = CordinatesRecover(ModelInfor, BestY, TaskInfor);  % get the cordinates of the best path among all the trials
        
        %% remark one : 需要修改
        FeasibleCostIndex=find(BestPathAndCost(1,:)<=3);
        FeasibleCost=BestPathAndCost(1, FeasibleCostIndex);
        SuccessfulRate=length(FeasibleCostIndex)/TrialTimes;
        MeanOfFeasible=mean(FeasibleCost);
        StdOfFeasible=std(FeasibleCost);
        
        
        MedianOfBestCost=median(BestPathAndCost(1,:));
        %MedianIndex = find(BestPathAndCost(1,:) == MedianOfBestCost);
        %MedianY=BestPathAndCost(2:Dimension+1, MedianIndex);
        %MedianBestPath=CordinatesRecover(ModelInfor, MedianY, TaskInfor);
        
        AverageConvergenceData=sum(ConvergenceData(:,2:end),1)/TrialTimes; % for the convergence graph
        
        fprintf('%i  %s  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f\r\n', TaskIndex, char(AlgorithmName(AlgorithmIndex)), SuccessfulRate, MinOfBestCost, MaxOfBestCost, MedianOfBestCost, MeanOfBestCost, StdOfBestCost, MeanOfFeasible, StdOfFeasible, AverageTime, StdTime);
        fprintf(find_file_result, '%i  %s  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f %.4f  %.4f  %.4f\r\n', TaskIndex, char(AlgorithmName(AlgorithmIndex)), SuccessfulRate, MinOfBestCost, MaxOfBestCost, MedianOfBestCost, MeanOfBestCost, StdOfBestCost, MeanOfFeasible, StdOfFeasible, AverageTime, StdTime);
        %fprintf(find_file_bestpath, '%i  %s  %.4f  %.4f\t\r\n',TaskIndex, char(AlgorithmName(AlgorithmIndex)), MinOfBestCost, BBestPath);
        %fprintf(find_file_bestpath, '\n%i  %s  %.4f\t', TaskIndex, char(AlgorithmName(AlgorithmIndex)), MinOfBestCost);
        fprintf(find_file_bestpath, '%.4f\t', BBestPath);
        
        fprintf(find_file_convergence_data, '%12.4e', AverageConvergenceData); %record the convergence data
        fprintf(find_file_eachTrialResult_data, '%12.4e', BestPathAndCost);
        fclose(find_file_convergence_data);
        fclose(find_file_eachTrialResult_data);
        fclose(find_file_bestpath);
    end
    
end
end
fclose ('all');
 