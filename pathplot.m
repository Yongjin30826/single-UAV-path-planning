function pathplot

Dimension=20;   % the number of waypoints
% path planning task information and parameters
% AlgorithmName={'PSO','CLPSOLSBFGS3','CLPSO'};
 AlgorithmName={'HFCLPSOLS'};
a=['-','*', '-p' ,'-+', '-d' ,'-s'];
for TaskIndex=1
    [~, ThreatInfor, ObstacleInfor ]=EnvironmentInfor(TaskIndex);
    [NT,~]=size(ThreatInfor);
    [NO,~]=size(ObstacleInfor);
    figure(TaskIndex);
    for AlgorithmIndex=[1]
        FileName=strcat(char(AlgorithmName(AlgorithmIndex)), 'Prob',int2str(TaskIndex), 'Dim',int2str(Dimension),'Path.txt');
        FindFile=fopen(FileName, 'r');
        Data=fscanf(FindFile,'%50f',[Dimension+2,inf]);
        X(AlgorithmIndex,:)=Data(:,1)';
        Y(AlgorithmIndex,:)=Data(:,2)';

        plot(X(AlgorithmIndex,:),Y(AlgorithmIndex,:),a(AlgorithmIndex*2-1:AlgorithmIndex*2));

        hold on
    end

    xlabel('x');
    ylabel('y');

    for ThreatIndex=1:NT
        r=ThreatInfor(ThreatIndex,3);
        xt=ThreatInfor(ThreatIndex,1);
        yt=ThreatInfor(ThreatIndex,2);
        theta = linspace(0,2*pi);
        x = r*cos(theta) + xt;
        y = r*sin(theta) + yt;
        plot(x,y,'b');
    end

    for ObstacleIndex=1:NO
        r=ObstacleInfor(ObstacleIndex,3);
        xt=ObstacleInfor(ObstacleIndex,1);
        yt=ObstacleInfor(ObstacleIndex,2);
        theta = linspace(0,2*pi);
        x = r*cos(theta) + xt;
        y = r*sin(theta) + yt;
        plot(x,y,'r');
    end
    
%     legend('PSO','CLPSOLSBFGS3','CLPSO');
legend('HFCLPSOLS');

%     close all;
    fclose ('all');
end
end
