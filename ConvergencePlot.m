function ConvergencePlot()

% AlgorithmName={'PSO','CLPSOLSBFGS3','CLPSO','RODDPSO', 'CLPSOLSBFGS','CLPSOLSBFGS2'};
AlgorithmName={'PSO','CLPSOLSBFGS3'};
a=['-*', '-p' ,'-+', '-d' ,'-s'];
for Dimension=20
    for ProblemIndex=[1]
        for AlgorithmIndex=1:2
            FileName=strcat(char(AlgorithmName(AlgorithmIndex)),'Prob',int2str(ProblemIndex),'Dim',int2str(Dimension),'Data.txt');
            FindFile=fopen(FileName, 'r');
            Result=fscanf(FindFile,'%50f',[1,inf]);
            Result=log10(Result);
            Length=length(Result);
            step=ceil(Length/10);
            figure(ProblemIndex)
            plot(Result,a(AlgorithmIndex*2-1:AlgorithmIndex*2),'MarkerIndices',1:step:Length);
            hold on;
        end
        xlabel('Function Evaluations');
        ylabel('log(cost)');
        %     legend('PSO','CLPSO','RODDPSO', 'CLPSOLS','CLPSOLLS');
        legend('PSO','CLPSOLSBFGS3');
        fclose all;
    end
end
end
