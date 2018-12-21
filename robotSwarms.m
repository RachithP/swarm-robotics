function robotSwarmFormation()
% Constants 
global N dT  KU_COEFF KS1_COEFF KS2_COEFF KT_COEFF

N = 5;         % number of bots
dT = .01;     % timestep length (position is updated each step)

% Control gains - This is taken from the 2004 paper. These are the gains
% for the individual abstract state variable used in the simulation of the
% paper
KU_COEFF = [2 0; 0 2];
KS2_COEFF = 2;
KT_COEFF = 2;

% Counter to accumulate the values of the data after each time interval
PLOT_COUNTER = 1;


%%%%%%%%%%%%%
% Main loop %
%%%%%%%%%%%%%
outpath = pwd;
% outputVideo = VideoWriter(fullfile(outpath,'SimulationVideo.mp4'),'MPEG-4');
% open(outputVideo);

% Place the robot swarm in the space 
bots = distributeBots(N);

% Time and posPlot to accumulate values
Time(1) = 0;
KS1_COEFF = 2; 
posPlot(N,1) = struct;

for i=1:N
    posPlot(i).qx(1) = bots(i).q(1) ;
    posPlot(i).qy(1) = bots(i).q(2) ;
    posPlot(i).uStarX(1) = 0 ;
    posPlot(i).uStarY(1) = 0 ;
    posPlot(i).uCapX(1) = bots(i).u(1) ;
    posPlot(i).uCapY(1) = bots(i).u(2) ;
end

PLOT_COUNTER = PLOT_COUNTER + 1;

%Initialize the abstract space 
uCentroid = [0 0].';
theta = 0;
s1 = 0;
s2 = 0;

%Initialize the matrices to calculate shape variables
E1 = [0 1; 1 0];
E2 = [1 0; 0 -1];
E3 = [0 -1; 1 0];


%Initializing the desired position variables 
uCentroidD = [3 ; 3];
s1D = 0.25;
s2D = 0.15;
thetaD = 0.9 ;

%Distance between bots - Taken from simulation of the current research
%paper
botRadius = 0.15;
botAxleLength = 0.1;
safeDist = 0.1 ;
sepDist = 2*(botRadius + botAxleLength) + safeDist;

%Initialize the configuration
[uCentroid,theta,s1,s2] = abstractSpace(bots);

figure;
drawElipseBoundary(bots,uCentroid,theta,s1,s2,uCentroidD,thetaD,s1D,s2D,botRadius, botAxleLength);
%Create a figure handle which we like to capture as a movie 
figure;
% F = getframe(gcf);

%Initial position of the robots with the elipse
drawElipseBoundary(bots,uCentroid,theta,s1,s2,uCentroidD,thetaD,s1D,s2D,botRadius, botAxleLength);
% writeVideo(outputVideo,getframe(gcf));

%Initialize plot variables
uPlotx(1) = uCentroidD(1) - uCentroid(1) ;
uPloty(1) = uCentroidD(1) - uCentroid(2) ;
thetaPlot(1) = thetaD - theta ;
s1Plot(1) = s1D - s1;
s2Plot(1) = s2D - s2;

% This variable represents the complete state of the system
xPlot(1) = norm([uPlotx(1) ; uPloty(1) ; thetaPlot(1) ; s1Plot(1) ; s2Plot(1)],2);

% step counter for every intervals of time
stepCounter = 0;
keepLooping = true;

%Values to contain the final velocity calculated from the convex
%optimization problem
uxMax = [0.1 ; 0.1];
kVel = [0 ; 0];
nVel = 0;

%Values to contain the min energy control law velocity. This is in
%accordance with the paper where the min energy control law velocity is
%having magnitude maximum of 0.05.
uxMaxCtrlLaw = 0.05;

conditionFailure = 0 ;

while (true == keepLooping && (350*1/dT >= stepCounter))
    
    % Calculating the abstract state variables
    [uCentroid,theta,s1,s2] = abstractSpace(bots);
    
    R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    H1 = eye(2) + R^2*E2;
    H2 = eye(2) - R^2*E2;
    H3 = R^2*E1 ;
    g = [[R uCentroid];0 0 1];
    
    %Check if the desired formation has been reached
    if(isequal(round((uCentroidD(1) - uCentroid(1)),3),0) && isequal(round((uCentroidD(2) - uCentroid(2)),3),0) && isequal(round((s1D-s1),3),0) && isequal(round((s2D-s2),3),0) && isequal(round((thetaD-theta),4),0))
        break;
    end
    
    %Calculate the error gains for each of the output
    dCentroid = KU_COEFF*(uCentroidD - uCentroid);
    dTheta = KT_COEFF*(thetaD - theta);
    dS1 = KS1_COEFF*(s1D - s1);
    dS2 = KS2_COEFF*(s2D - s2);
    
    for robotInd1=1:N
        position = bots(robotInd1).q.'; % Current position of the robot
        %Calculation of velocity using min energy control law
        velocity = dCentroid + ((s1-s2)*H3*(position - uCentroid)*dTheta / (s1+s2))....
            + (H1*(position - uCentroid)*dS1 / 4*s1) + (H2*(position - uCentroid)*dS2/4*s2);
        % ui*
        
        %Scaling the values of the min energy ctrl velocity to 0.05ms-1
        nVelCtrlLaw = max(1,norm(velocity,2)/uxMaxCtrlLaw);
        
        velocity = velocity / nVelCtrlLaw ;
        
        % Converting u to v using R and also position w.r.t moving frame
        movPosition = R.'*(position - uCentroid); %pi
        movVelocity = R.'*velocity; %vi*
        
        %Inequality constraint for asymptotic convergence
        % stateTilde is the error of the state
        stateTilde = [uCentroidD - uCentroid;thetaD-theta;s1D-s1;s2D-s2];
        
        % Gamma is the transformation matrix from moving frame to abstract
        % space
        Gamma = [g zeros(3,2);zeros(2,3) eye(2)];
        
        % 5 by 1 matrix used in the monotonic convergence criterian
        val = [eye(2);(1/s1-s2)*movPosition.'*E1;movPosition.'*(eye(2)+E2);movPosition.'*(eye(2)-E2)];
        
        %Gain Matrix - 5*5 matrix
        GainMat = [KU_COEFF(1,:) 0 0 0; KU_COEFF(2,:) 0 0 0; 0 0 KT_COEFF 0 0; 0 0 0 KS1_COEFF 0; 0 0 0 0 KS2_COEFF];
        
        Acondition1 = stateTilde.'*GainMat*Gamma*val;
        
        %Inequality constraint to saturate the maximum velocity of the robot
        % calculated using convex optimization to max of 0.1ms-1
        % u = Rv
        Acondition2 = [1 0] * R ;
        Acondition3 = [0 1] * R;
        
        AMatCondition =[-Acondition1 ; Acondition2 ; Acondition3];
        BMatCondition = [0; 0.1; 0.1];
        
        %Check for conditions when the robots are within collision distance.
        %This is when the ccollision avoidance constrained is applied for the
        %robots
        for robotInd2=1:N
            % calculate the position of each robot and compare against
            % current
            if(robotInd2 ~= robotInd1)
                bot1Position = R.'*(bots(robotInd1).qN.' - uCentroid); %p1
                bot2Position = R.'*(bots(robotInd2).qN.' - uCentroid); %p2 - Old
                                
                bot2Velocity = R.'*(bots(robotInd2).uN.');%v2 - New
                %the calculated delta value
                delta = norm(bot2Position - bot1Position, 2);
                
                if(delta <= sepDist)
                    % Condition for the collision avoidance
                    Acondition4 = (bot1Position - bot2Position).';
                    Bcondition4 = ((bot1Position - bot2Position).'*bot2Velocity);
                    AMatCondition = [AMatCondition ; -Acondition4];
                    BMatCondition =  [BMatCondition ; -Bcondition4];
                end
            end
        end
        
        % Solving for optimal velocity based on previous condition
        opts1 = optimset('display','off');
        velocityCap = lsqlin(sqrt(2)*eye(2),sqrt(2)*movVelocity,AMatCondition,BMatCondition,[],[],[],[],[],opts1); %vicap
        
        % There are instances in which the lsqlin function fails. This is a
        % check for the failure to debug the system
        if(round(AMatCondition*velocityCap,3) > round(BMatCondition,3))
            AMatCondition*velocityCap - BMatCondition
            conditionFailure = conditionFailure + 1;
        end
        
        % velocity with respect to the world frame
        velWorldFrame = R*velocityCap ; %uicap
        
        %As per the paper we will contain the velocity at max cap velocity
        %for the robot
        kVel(1) = max(1, abs(velWorldFrame(1))/uxMax(1));
        kVel(2) = max(1, abs(velWorldFrame(2))/uxMax(2));
        
        nVel = max(kVel(1),kVel(2));
        
        velWorldFrame = velWorldFrame / nVel ;
        
        % Update the position of the robot based on Euler method for simulation
        % Modelling the position based on the differential drive robot
        % model
        updateDifferentialDrivePos(velWorldFrame, botAxleLength, robotInd1);
        
        bots(robotInd1).uStar = velocity.' ; %ui*
    end
    
    %Update each of the robot position
    for robotInd1=1:N
        
        %Update the position and velocity variables of the robots
        bots(robotInd1).q = bots(robotInd1).qN ;
        bots(robotInd1).u = bots(robotInd1).uN;
        bots(robotInd1).tr = bots(robotInd1).trN ;
        
        % Accumulate the plot variables
        % Position of robots
        posPlot(robotInd1).qx(PLOT_COUNTER) = bots(robotInd1).q(1);
        posPlot(robotInd1).qy(PLOT_COUNTER) = bots(robotInd1).q(2);
        % Minimum energy control velocity
        posPlot(robotInd1).uStarX(PLOT_COUNTER) = bots(robotInd1).uStar(1);
        posPlot(robotInd1).uStarY(PLOT_COUNTER) = bots(robotInd1).uStar(2);
        % Control velocity input based on convex optimization
        posPlot(robotInd1).uCapX(PLOT_COUNTER) = bots(robotInd1).u(1);
        posPlot(robotInd1).uCapY(PLOT_COUNTER) = bots(robotInd1).u(2);
        
    end
    
    stepCounter = stepCounter+1;
    
    %Accumulate the data values for the plot
    Time(PLOT_COUNTER) = (Time(1) + dT*stepCounter);
    uPlotx(PLOT_COUNTER) = (uCentroidD(1) - uCentroid(1)) ;
    uPloty(PLOT_COUNTER) = (uCentroidD(2) - uCentroid(2)) ;
    thetaPlot(PLOT_COUNTER) = (thetaD - theta) ;
    s1Plot(PLOT_COUNTER) = (s1D-s1);
    s2Plot(PLOT_COUNTER) = (s2D - s2);
    xPlot(PLOT_COUNTER) = norm([uPlotx(PLOT_COUNTER) ; uPloty(PLOT_COUNTER) ; thetaPlot(PLOT_COUNTER) ; s1Plot(PLOT_COUNTER) ; s2Plot(PLOT_COUNTER)],2);
    
    PLOT_COUNTER = PLOT_COUNTER+1;
    
    if(0 == isnan(s1) && 0 == isnan(s2) )
        s1 , s2, uCentroid, theta
    end
    
    if(mod(stepCounter,20) == 0)
        drawElipseBoundary(bots,uCentroid,theta,s1,s2,uCentroidD,thetaD,s1D,s2D,botRadius, botAxleLength); 
%         writeVideo(outputVideo,getframe(gcf));
    end

end

%Generate the data for difference in bot position
l=1;
for j=1:N
    for k=(j+1):N
        for plotCount=1:PLOT_COUNTER-1
            tempMatVal = [posPlot(j).qx(plotCount) - posPlot(k).qx(plotCount) ; posPlot(j).qy(plotCount) - posPlot(k).qy(plotCount)];
            posDiffPlot(l).q(plotCount) = norm(tempMatVal,2);
            
        end
        l = l+1;
    end
end

%Generate the data for the velocity plot
for j=1:N
    for plotCount=1:PLOT_COUNTER-1
    velStarPlot(j).u(plotCount) = norm([posPlot(j).uStarX(plotCount) posPlot(j).uStarY(plotCount)],2);
    velCapPlot(j).u(plotCount) = norm([posPlot(j).uCapX(plotCount) posPlot(j).uCapY(plotCount)],2);
    end
end

% Draw the final elipse position
drawElipseBoundary(bots,uCentroid,theta,s1,s2,uCentroidD,thetaD,s1D,s2D,botRadius, botAxleLength);
% writeVideo(outputVideo,getframe(gcf));

% close(ouVtputVideo);

%We start Plotting our parameters here
% Plot of the error in the state variables
figure ;

plot1 = plot(Time,xPlot,Time, abs(uPlotx),'--',Time, abs(uPloty),'--',Time, abs(thetaPlot),'--', Time, abs(s1Plot), '--',Time, abs(s2Plot), '--','linewidth', 0.7);
plot1(1).LineWidth = 1;
plot1(1).Color = 'black';
xlabel('Time(s)');
ylabel('Magnitude');
legend('||x||','||x_1||','||x_2||','||x_3||','||x_4||','||x_5||');
title('convergencePlot');
hold on;


figure ;
plot(Time,2*(botRadius + botAxleLength)*ones(size(Time)),'--','Color','black');
hold on;
for i=1:l-1
    plot(Time,posDiffPlot(i).q);
    hold on;
end
xlabel('Time(s)');
ylabel('||q_i-q_j||(m)');
title('RobotPositionComparison');

hold on;

figure ;
plot(Time, velStarPlot(1).u,'--', Time, velCapPlot(1).u,'linewidth',1);
xlabel('Time(s)');
ylabel('Magnitude(m/s)');
title('VelocityPlot');
legend('||u_i^{*}||','||u_i^{\^}||');
hold on;

figure;
drawElipseBoundary(bots,uCentroid,theta,s1,s2,uCentroidD,thetaD,s1D,s2D,botRadius,botAxleLength);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bot initialization                                                      % 
% Description - This function distributes the robots in Euclidean space as%
% a normal distribution. The separation between each robot should be      %
% greater than the initial configuration                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function bots = distributeBots(N)
        
        %Structure variables: 
        %q - holds the position w.r.t the world frame
        %u - Velocity w.r.t the world frame
        %qN - Holds the new position of the robots that has to be updated
        bots(N,1) = struct('q',[0; 0],'u',[0; 0].','qN',[0 ; 0],'uN',[0 ; 0],'uStar',[0 ; 0],'tr',0,'trN',0);
        
        
        % Choosing Random distribution of the robots where the robots are
        % separated at a distance greater thatn the safe separation distance
        % between each of them. We choose
        % an arbitrary value of the mean and the standard deviation
        
        randLoop = true;
        while(randLoop)
            counter = 0;
            % Standard deviation is equal to robot count / 4
            A = normrnd(5,max(1,N/4),[2,N]);
            
            for robotInd1=1:N
                for robotInd2=1:N
                    if(robotInd1 ~= robotInd2)
                        if(norm((A(:,robotInd1) - A(:,robotInd2)),2) < 0.6)
                            counter = counter+1;
                        end
                    end
                end
            end
            
            if(counter == 0)
                %Print the position of the robots
                A
                randLoop = false;
            end
        end
        
        for robotInd=1:N
            % place agent, Initializing all states of the robots
            bots(robotInd).q = A(1:2,robotInd).';
            bots(robotInd).qN = A(1:2,robotInd).';
            bots(robotInd).u = [0 0];
            bots(robotInd).uN = [0 0];
            bots(robotInd).uStar = [0 0];
            bots(robotInd).tr = 0 ;
            bots(robotInd).trN = 0 ;
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name - Abstract Space                                                   % 
% Description - Computes the abstract state variables of the formation    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [aCentroid,aTheta,aS1,aS2] = abstractSpace(bots)
        E1 = [0 1; 1 0];
        E2 = [1 0; 0 -1];
        E3 = [0 -1; 1 0];
        
        tThetaY = 0;
        tThetaX = 0;
        tS1 = 0;
        tS2 = 0;
        tempCentroid = 0;
        
        %Calculate Centroid
        for robotInd =1:N
            position = bots(robotInd).q.';
            tempCentroid = tempCentroid + position ;
        end
        aCentroid = 1/N * tempCentroid ;
        
        %Calculate Orientation
        for robotInd = 1:N
            position = bots(robotInd).q.';
            tThetaY = tThetaY + (position - aCentroid).'*E1*(position - aCentroid);
            tThetaX = tThetaX + (position - aCentroid).'*E2*(position - aCentroid);
        end
            aTheta = atan2(tThetaY , tThetaX)/2 ;
            R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
            H1 = eye(2) + R^2*E2;
            H2 = eye(2) - R^2*E2;
            H3 = R^2*E1 ;
        
        for robotInd=1:N
            position = bots(robotInd).q.';
            tS1 = tS1 + ((position - aCentroid).'*H1*(position - aCentroid)) ;
            tS2 = tS2 + ((position - aCentroid).'*H2*(position - aCentroid)) ;
        end
        
        aS1 = tS1 / (2*(N-1)) ;
        aS2 = tS2 / (2*(N-1)) ;

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name - drawElipseBoundary                                               % 
% Description - Draw elipse boundary around the robot formation and the   % 
% final desired positon                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function drawElipseBoundary(aBots, aCentroid , aTheta, aS1, aS2 , aCentroidD, aThetaD, aS1D, aS2D, botRadius,botAxleLength)
        % Plot the robot positions and the ellipsoid
        X = [];
        t=0:0.01:2*pi;
        
        for robotInd=1:N
            X = [X [aBots(robotInd).q(1);aBots(robotInd).q(2)]];
        end
        
        if aS1>aS2
            % Multiply [acos(t); bsin(t)] by R. Taking confidence parameter as 1
            x1 = aCentroid(1) + sqrt(9.2103*aS1)*cos(t)*cos(aTheta) - sqrt(2*aS2)*sin(t)*sin(aTheta);
            y1 = aCentroid(2) + sqrt(9.2103*aS2)*sin(t)*cos(aTheta) + sqrt(2*aS1)*cos(t)*sin(aTheta);
        else
            % Multiply [acos(t); bsin(t)] by R. Taking confidence parameter as 1
            x1 = aCentroid(1) + sqrt(9.2103*aS2)*cos(t)*cos(aTheta) - sqrt(2*aS1)*sin(t)*sin(aTheta);
            y1 = aCentroid(2) + sqrt(9.2103*aS1)*sin(t)*cos(aTheta) + sqrt(2*aS2)*cos(t)*sin(aTheta);
        end
        
        %Plot the desired elipse position
                
        if aS1D>aS2D
            % Multiply [acos(t); bsin(t)] by R. Taking confidence parameter as 1
            x2 = aCentroidD(1) + sqrt(9.2103*aS1D)*cos(t)*cos(aThetaD) - sqrt(2*aS2D)*sin(t)*sin(aThetaD);
            y2 = aCentroidD(2) + sqrt(9.2103*aS2D)*sin(t)*cos(aThetaD) + sqrt(2*aS1D)*cos(t)*sin(aThetaD);
        else
            % Multiply [acos(t); bsin(t)] by R. Taking confidence parameter as 1
            x2 = aCentroidD(1) + sqrt(9.2103*aS2D)*cos(t)*cos(aThetaD) - sqrt(2*aS1D)*sin(t)*sin(aThetaD);
            y2 = aCentroidD(2) + sqrt(9.2103*aS1D)*sin(t)*cos(aThetaD) + sqrt(2*aS2D)*cos(t)*sin(aThetaD);
        end
        
%         labels = {'1','2','3','4','5','6','7','8','9','10'};
        plot(X(1,:),X(2,:),'o');
%         text(X(1,:),X(2,:),labels,'VerticalAlignment','bottom','HorizontalAlignment','right');
        title('Robot Configuration')
        xlabel('X-axis')
        ylabel('Y-axis')
        hold on
        plot(x1,y1,'r',x2,y2,'b');
        for robotInd=1:N                    
            botCenter = aBots(robotInd).q;
            botOrientation = aBots(robotInd).tr;
            rotMat = [cos(botOrientation) -sin(botOrientation);sin(botOrientation) cos(botOrientation)];
%             rotMat = [cos(botOrientation) sin(botOrientation) botRadius*cos(botOrientation); -sin(botOrientation) cos(botOrientation) botRadius*sin(botOrientation); 0 0 1];
            lineEnd = botCenter + [(botRadius+botAxleLength)*cos(botOrientation) (botRadius+botAxleLength)*sin(botOrientation)];
%             lineEnd = lineEnd(1:2,1).';
            % Defining circles around each robot
            x = aBots(robotInd).q(1)+ botAxleLength*cos(botOrientation) + (botRadius+botAxleLength)*cos(t);
            y = aBots(robotInd).q(2)+ botAxleLength*sin(botOrientation) + (botRadius+botAxleLength)*sin(t);
            plot(x,y,'g');
            quiver(botCenter(1,1),botCenter(1,2),lineEnd(1,1)+botAxleLength*cos(botOrientation) - botCenter(1,1), lineEnd(1,2)+botAxleLength*sin(botOrientation) - botCenter(1,2),0,'Color','red');           
            
        end
        legend('Robot Position','Actual Formation','Desired Formation');
        hold off
    end

    function updateDifferentialDrivePos(velWorldFrame, axleLength, robotInd)
        linVel = 0 ;
        angVel = 0 ;
        
        botOrientation = bots(robotInd).tr ; %New Orientation
        rotMat = [cos(botOrientation) sin(botOrientation) ; -sin(botOrientation)/axleLength cos(botOrientation)/axleLength];
        
        velMat = rotMat*velWorldFrame ;
        
        linVel = velMat(1);
        angVel = velMat(2);
        
        %Update the kinematic position of the robots. We use the same model
        %that is used in the robotic simulator toolbox to update the
        %position
        dx = dT * linVel * cos(botOrientation);
        dy = dT * linVel * sin(botOrientation);
        
        dtr = dT * angVel ;
        
        %Update the new positions and new velocities to which the robots
        %have to move and the new orientation of the robots
        bots(robotInd).qN(1)  = bots(robotInd).q(1) + dx ;
        bots(robotInd).qN(2) = bots(robotInd).q(2) + dy;
        
        bots(robotInd).uN = velWorldFrame.'; 
        
        bots(robotInd).trN = bots(robotInd).tr + dtr ;
    end

end