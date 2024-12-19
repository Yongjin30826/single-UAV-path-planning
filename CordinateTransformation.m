function  Infor= CordinateTransformation(TaskInfor, ThreatInfor,ObstacleInfor, Num_WayPoints)
      StartPoint=TaskInfor(1:2);
      TargetPoint=TaskInfor(3:4);
      %UpperBound=[0, TaskInfor(5)];
      %LowerBound=[TaskInfor(6), 0];
      d=dist(StartPoint, TargetPoint');  % 欧式距离公式，二维两点之间的距离
      Theta=asin((TargetPoint(2)-StartPoint(2))/d);  % 对应论文中的原直角坐标系的x轴逆时针旋转到新坐标轴的角度
      a=[cos(Theta) sin(Theta); -sin(Theta), cos(Theta)];
      % transform task cordinates
      TStartPoint=a*(StartPoint-StartPoint)';
      TTargetPoint=a*(TargetPoint-StartPoint)'; % (2,1)
      TTaskInfor(1:2)=TStartPoint';  % 矩阵转置
      TTaskInfor(3:4)=TTargetPoint';
      TTaskInfor(5:6)=TaskInfor(5:6); 
      % transform threat cordinates   %  转换威胁点的坐标
      [n, ~]=size(ThreatInfor);  % 返回值n为threat的行，即代表了威胁点的数量
      TThreatInfor=ThreatInfor;
      for i=1:n
          TThreatInfor(i,1:2)=(a*(ThreatInfor(i, 1:2)-StartPoint)')';
      end
      
      k=1:Num_WayPoints;       
      x=k.*(d/(Num_WayPoints+1)); % (1, Num_WayPoints)  D的数量，将线段分成D+1段
      % transform obstable cordinates
      [m, ~]=size(ObstacleInfor);   % 返回值m代表obstacleinfor的行，即代表了障碍物的数量
      TObstacleInfor=ObstacleInfor;  
      for i=1:m
          TObstacleInfor(i,1:2)=(a*(ObstacleInfor(i, 1:2)-StartPoint)')';
      end
      
      % calculate the transformed bound   计算换算的边界值
      kk=(TargetPoint(2)-StartPoint(2))/(TargetPoint(1)-StartPoint(1));  % 原坐标ST直线的斜率
      % y=k*(x-StartPoint(1))+StartPoint(2)
      for j=1:Num_WayPoints
          xx(j)=StartPoint(1)+cos(Theta)*x(j); % the x  in the original cordinate   点在原直角坐标系下的横坐标= 新的直角坐标系原点的横坐标 + 点在新的直角坐标下的坐标✖余弦值
          yy(j)=kk*(xx(j)-StartPoint(1))+StartPoint(2);  % 点在原直角坐标系下的纵坐标
          b(j)=yy(j)+(1/kk)*xx(j);     % the intersection with y axis  两个y轴的交点
          if b(j)<=TaskInfor(6)
              UpperBound(j)=b(j);  % 上界
              UpperPoint(j,:)=[0, UpperBound(j)];
          else
              UpperPointX(j)=-(TaskInfor(6)-b(j))*kk;
              UpperPoint(j,:)=[UpperPointX(j), TaskInfor(6)];
          end
          TUpperPoint(j,:)=a*(UpperPoint(j,:)-StartPoint)';
          
          bb(j)= kk*b(j);       % the intersection with x axis  两个x轴的交点
          if bb(j)<=TaskInfor(5)
              LowerBound(j)=bb(j);
              LowerPoint(j,:)=[LowerBound(j), 0];
          else
              LowerPointY(j)=-(1/kk)*TaskInfor(5)+b(j);
              LowerPoint(j,:)=[TaskInfor(5), LowerPointY(j)];
          end
          TLowerPoint(j,:)=a*(LowerPoint(j,:)-StartPoint)';
      end
    
      
      
      %b=a*(UpperBound-StartPoint)';
      %c=a*(LowerBound-StartPoint)';
      Bound=[TLowerPoint(:,2),TUpperPoint(:,2)]; % lowerBound, UpperBound  (D,2) 只是针对新的坐标系下的y轴的范围
      %Bound=[-100, 100];
      Infor.Task=TTaskInfor;
      Infor.Threat=TThreatInfor;
      Infor.Obstacle=TObstacleInfor;
      Infor.x=x';
      Infor.Bound=Bound;
      Infor.Num_WayPoints=Num_WayPoints;
      Infor.a=a;
end