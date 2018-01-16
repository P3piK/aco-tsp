clear
clc

%cities 2
x = [0 2 6 7 15 12 14 9.5 7.5 0.5];
y = [1 3 5 2.5 -0.5 3.5 10 7.5 9 10];

%cities 3
%x = [0 3 6 7  15   10 16 5 8 1.5];
%y = [1 2 1 4.5 -1 2.5 11 6 9 12];

% Initialization
m = 100; % ants
N = 10; % cities
numberOfIterations = 100;

DistanceTraveled = Inf([numberOfIterations m]);

alpha = 1;
beta = 5;
rho = 0.5;

M = zeros([1 N]); % working memory
a = zeros([ N N numberOfIterations]);
Tau = zeros([N N numberOfIterations]);
T0 = 0.3;
DeltaTau = zeros([N N]);
Tau(:,:,1) = (1-rho) * T0;

% distance between cities
for i = 1:(N-1)
    for j=i:N
        Distance(i, j) = (x(i)-x(j))^2 + (y(i)-y(j))^2;
        % Distance(i, j) = randi(100);
    end
end
for i = 2:N
    for j = 1:(i - 1)
        Distance(i,j) = Distance(j,i);
    end
end

% initializing Probability table + zeroing diagonal in Distance
for i = 1:N
    for j = 1:N
        Probability(i,j) = 1 / (N-1);
    end
    Probability(i,i) = 0;
end

% main loop
for t=1:numberOfIterations
    for k=1:m
        DistanceTraveled(t,k) = 0;
        startCity = randi(N);
        M(1) = startCity;
        tempProbability = Probability;
        
        % excluding this city from available cities
        for zeroProbability=1:N
            tempProbability(zeroProbability, M(1)) = 0;
        end
        
        % choosing path
        for n=2:N
                      
            % generate random to choose
            nextCity = 0;
            while nextCity < 1 || nextCity > N
                r = rand;
                nextCity = sum(r >= cumsum([0, tempProbability(M(n-1),:)])); % choose with probability
            end

            
            M(n) = nextCity;  % add city to working memory
            DistanceTraveled(t, k) = DistanceTraveled(t, k) + Distance(M(n-1), M(n));
   
            % decision table A(k)=a(ij)(t)        
            % denominator
            tempDenominatorSum = 0;
            for tempDenominatorVar1 = 1:N
                if tempProbability(M(n-1),tempDenominatorVar1) ~= 0
                    tempDenominatorSum = tempDenominatorSum + (Tau(M(n-1), tempDenominatorVar1,t))^alpha * (1 / Distance(M(n-1),tempDenominatorVar1))^beta;
                end
            end
            % num / denom
            a(M(n-1), M(n), t) = a(M(n-1),M(n),t) + (Tau(M(n-1),M(n),t))^alpha * (1/Distance(M(n-1),M(n)))^beta / tempDenominatorSum;
            
            % exclude this city from available cities
            for zeroProbability=1:N
                tempProbability(zeroProbability, M(n)) = 0;
            end
                        
            % update probabilities
            probabilityCount = 0;
            for cityCount=1:N
                if tempProbability(nextCity,cityCount) ~= 0
                    probabilityCount = probabilityCount + tempProbability(nextCity, cityCount);
                end
            end
            for nonzeroProbabilityCity=1:N
                if tempProbability(nextCity, nonzeroProbabilityCity) ~= 0
                    tempProbability(nextCity, nonzeroProbabilityCity) = tempProbability(nextCity, nonzeroProbabilityCity) / probabilityCount;
                end
            end
        end
        
        
        DistanceTraveled(t,k) = DistanceTraveled(t,k) + Distance(M(N), M(1)); % route length
        % quantity of pheromone on the travelled arc
        for tempDeltaTauVar = 2:N
            firstCity = M(tempDeltaTauVar - 1);
            secondCity = M(tempDeltaTauVar);
            DeltaTau(firstCity, secondCity) = DeltaTau(firstCity, secondCity) + (1 / DistanceTraveled(t,k));
        end
        
        DeltaTau(M(N), M(1)) = DeltaTau(M(N), M(1)) + (1 / DistanceTraveled(t,k));
        Tau(:,:,t+1) = (1-rho) * Tau(:,:,t) + DeltaTau;
        
        % plot shortest path
        if DistanceTraveled(t,k) == min(DistanceTraveled(:))
            clf
            for drawVar=1:N
               rectangle('Position', [x(M(drawVar))-0.1, y(M(drawVar))-0.1, 0.2, 0.2], 'FaceColor', [0, .5, .5])
               if drawVar~=N
                   line([x(M(drawVar)) x(M(drawVar+1))],[y(M(drawVar)) y(M(drawVar+1))])
                   hold on
               else
                   line([x(M(drawVar)) x(M(1))], [y(M(drawVar)) y(M(1))])
                   
               end
            end
            drawnow
        end
        
    end
    
    % update probability from decision table
    for Row=1:N
        for Col=1:N
            Probability(Row,Col) = a(Row,Col,t) / sum(a(Row,:,t));
            if Probability(Row,Col) == 0
                Probability(Row,Col) = 0.00001;
            end
        end
    end
    
    
end

Distance
Probability
minDist = min(DistanceTraveled(:))
% [row,col] = find(DistanceTraveled==minDist)
