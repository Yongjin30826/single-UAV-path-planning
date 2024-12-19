function FitValue=CostFunction(ParticleSwarm, ModelInfor)

[Dimension, SwarmSize]=size(ParticleSwarm); % Dimension is the number of waypoints
X=ModelInfor.x;  % x是新的直角坐标系下的坐标
Threat=ModelInfor.Threat;
Obstacle=ModelInfor.Obstacle;
Task=ModelInfor.Task;
Penalty=3;  % 论文中的碰撞P
if Task(5)<=100 
    MaximumLength=300;    % 不懂
else
    MaximumLength=500;     % 不懂
end
    
StartPoint=Task(1:2);
TargetPoint=Task(3:4);
ST=dist(StartPoint, TargetPoint');  % the length of the straight line connencting the starting point and the target point

for i=1: SwarmSize
   Y=ParticleSwarm(:, i);
   if ~isreal (Y)
       Y=real(Y);
       %FitValue(i)=nan;
       %break;
   end

   Waypoints=[X Y];
   Path=[StartPoint; X Y; TargetPoint]; % (D+2,2)
   
   %Calculate the cost associated with the total length \in [0,1]
   %长度成本，1-Lsg/Lpath
   d=Path(2:Dimension+2,:)-Path(1:Dimension+1,:);  % A（i:i+m,:）， 第i~i+m行的全部元素  d=[x增量，y增量]
   PathLength=sum((sum(d.*d,2)).^0.5); %sum(x,2)对每一行分别求和
   CLength=1-(ST/PathLength);
   
   %Another way to calculate the cost associated with the threats \in [0,1]
   %计算威胁成本
if ~isempty (Threat)
   [m, ~]=size(Threat);  % m: the number of threats  威胁点的数量
   dThreat=pdist2(Path, Threat(:,1:2)); % dThreat(i,j) is the distance between i (D+2) point and j threat 
   % 计算i(D+2)个点与j威胁点之间的距离
  
   for j=1:m  % for each threat
       for k=2:Dimension+2  % for each point/line segment
          if dThreat(k,j)<=Threat(j,3)   
              % j是威胁点的编号数量，矩阵threat的第三列则是威胁点的半径。此为编号从2到D+2个点到威胁点的距离小于威胁半径
               if dThreat(k-1,j)<=Threat(j,3)  % the segment is in the circle
                   % 如果前一个点也在威胁圈里面，说明两点之间会在这个威胁圈里
                
                   TS(k)=norm(Path(k,:)-Path(k-1,:));    % 线段在威胁圈里边
                   
                 % norm的向量运算，计算威胁圈里两个点之间的欧氏距离
               else                            % intersect 相交
                   x1=Path(k-1,1); y1=Path(k-1,2);
                   x2=Path(k,1); y2=Path(k,2);      
                   x3=Obstacle(j,1);y3=Obstacle(j,2);r=Obstacle(j,3);
                   A=(x2-x1)^2+(y2-y1)^2;
                   B=2*((x2-x1)*(x1-x3)+(y2-y1)*(y1-y3));
                   C=x3^2+y3^2+x1^2+y1^2-2*(x3*x1+y3*y1)-r^2;
                   delta=B^2-4*A*C;
                   mu1=(-B+delta^0.5)/(2*A); mu2=(-B-delta^0.5)/(2*A);
                   if  0<=mu1&& mu1<=1
                       
                       TS(k)=norm(Path(k,:)-Path(k-1,:))*(1-mu1);   %前一个点在威胁圈外，后一个点在威胁圈内
                   else
                       if  0<=mu2 && mu2<=1
                           
                           TS(k)=norm(Path(k,:)-Path(k-1,:))*(1-mu2);    %前一个点在威胁圈外，后一个点在威胁圈内
                       else
                           TS(k)=0;                                      %线段一个点在威胁圈外，后一个点在威胁圈上（边）
                       end
                   end
                   
               end
          else 
              if dThreat(k-1,j)<=Threat(j,3)   % intersect
                  x1=Path(k-1,1); y1=Path(k-1,2);
                  x2=Path(k,1); y2=Path(k,2);
                  x3=Obstacle(j,1);y3=Obstacle(j,2);r=Obstacle(j,3);
                  A=(x2-x1)^2+(y2-y1)^2;
                  B=2*((x2-x1)*(x1-x3)+(y2-y1)*(y1-y3));
                  C=x3^2+y3^2+x1^2+y1^2-2*(x3*x1+y3*y1)-r^2;
                  delta=B^2-4*A*C;
                  mu1=(-B+delta^0.5)/(2*A); mu2=(-B-delta^0.5)/(2*A);
                  if  0<=mu1&& mu1<=1
                       TS(k)=norm(Path(k,:)-Path(k-1,:))*mu1;    % 前一个点在威胁圈内，后一个点在威胁圈外
                   else
                       if  0<=mu2 && mu2<=1
                           TS(k)=norm(Path(k,:)-Path(k-1,:))*mu2;  %前一个点在威胁圈内，后一个在威胁圈外
                       else
                           TS(k)=0;     % 前一个点在威胁圈上（边），后一个点威胁圈外
                       end
                   end
              else
                  x1=Path(k-1,1); y1=Path(k-1,2);
                  x2=Path(k,1); y2=Path(k,2);
                  x3=Obstacle(j,1);y3=Obstacle(j,2);r=Obstacle(j,3);
                  A=(x2-x1)^2+(y2-y1)^2;
                  B=2*((x2-x1)*(x1-x3)+(y2-y1)*(y1-y3));
                  C=x3^2+y3^2+x1^2+y1^2-2*(x3*x1+y3*y1)-r^2;
                  delta=B^2-4*A*C;
                  if delta<=0      % no intersection
                      TS(k)=0;
                  else             % two intersections 
                      mu1=(-B+delta^0.5)/(2*A); mu2=(-B-delta^0.5)/(2*A);
                      
                     TS(k)=norm(Path(k,:)-Path(k-1,:))*(mu1-mu2);            %前一个点在威胁圈外，后一个点也在威胁圈外
                     
                  end
              end
          end
       end
   end
   CDanger=sum(TS)/sum(Threat(:,3));
   if CDanger>1
       CDanger=1;
   end
else
    CDanger=0;
end
   
   %Calculate the cost associated with turning \in [0,1]   计算转弯成本
   for ii=2:Dimension+1
       %Theta=d(ii,:)*d(ii-1,:)'/(d(ii,:)*d(ii,:)'^0.5)*((d(ii,:)*d(ii,:)'^0.5));
       turning(ii-1)=dot(d(ii,:),d(ii-1,:))/(norm(d(ii,:))*norm(d(ii-1,:)));  % dot实数向量的点积
   end
   CTurning=(1-mean(turning))/2;  % 范围也在0-1之间
   
   
   
   %Calculate the cost associated with the collision with obstacles \in [p,p+1]
   %计算碰撞成本
   [n, ~]=size(Obstacle);  % n: the number of obstacles 
   dObstacle=pdist2(Path, Obstacle(:,1:2)); %dd(i,j) is the distance between i (D+2) point and j obstacle 
   % path上的点到障碍物中心点的距离
   Collision=0;
   for j=1:n  % for each obstacle
       for k=2:Dimension+2  % for each point. Note: the start point is collision free
          if dObstacle(k,j)<Obstacle(j,3)
               Collision=Collision+1;
          else 
          x1=Path(k-1,1); y1=Path(k-1,2);
          x2=Path(k,1); y2=Path(k,2);
          x3=Obstacle(j,1);y3=Obstacle(j,2);r=Obstacle(j,3);
          A=(x2-x1)^2+(y2-y1)^2;
          B=2*((x2-x1)*(x1-x3)+(y2-y1)*(y1-y3));
          C=x3^2+y3^2+x1^2+y1^2-2*(x3*x1+y3*y1)-r^2;
          delta=B^2-4*A*C;
          if delta<=0
              % do nothing
          else
              mu1=(-B+delta^0.5)/(2*A); mu2=(-B-delta^0.5)/(2*A);
              if (mu1<1 && mu1>0 ) || (mu2<1 && mu2>0)
                  Collision=Collision+1;
              end
          end
          end
       end
   end
   if Collision>0
       CCollision=Penalty+Collision/(Dimension+1);
   else
       CCollision=0;
   end
   
   %Calculate the cost associated with the fuel represented by fly length \in [p,p+1]
   if PathLength<=MaximumLength
       CFuel=0;
   else
       CFuel=Penalty+(PathLength-MaximumLength)/MaximumLength;
   end
   
   %Calculate the cost associated with the smoothing \in [p,p+1]

  %CLength
  %CDanger
  %CTurning
  %CCollision
  %CFuel
   FitValue(i)=CLength+CDanger+CTurning+ CCollision+CFuel;
   
end

end