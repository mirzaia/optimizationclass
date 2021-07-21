% Kelas Optim 4 Maret 2020
% Tugas: Simulasi, lalu di modifikasi


% Problem: Travelling Salesman
% 1. Draw the Maps and Stops
figure;

load('usborder.mat','x','y','xx','yy');
rng(10,'multFibonacci') % makes a plot with stops in Maine & Florida, and is reproducible METODE Multiplicative Lagged Fibonacci'
nStops = 69; % you can use any number, but the problem size scales as N^2 JUMLAH PEMBERHENTIAN
stopsLon = zeros(nStops,1); % allocate x-coordinates of nStops
stopsLat = stopsLon; % allocate y-coordinates
n = 4; %diubah
while (n <= nStops)
    xp = rand*3; %diubah
    yp = rand;
    if inpolygon(xp,yp,xx,yy) % test if inside the border
        stopsLon(n) = xp;
        stopsLat(n) = yp;
        n = n+1;
    end
end
plot(x,y,'Color','yellow'); % draw the outside border WARNA BORDER 
hold on
% Add the stops to the map
plot(stopsLon,stopsLat,'*b')
hold off

% Problem Formulation

% Formulate the travelling salesman problem for integer linear programming as follows:

% Generate all possible trips, meaning all distinct pairs of stops.
% Calculate the distance for each trip.
% The cost function to minimize is the sum of the trip distances for each trip in the tour.
% The decision variables are binary, and associated with each trip, where each 1 represents a trip that exists on the tour, and each 0 represents a trip that is not on the tour.
% To ensure that the tour includes every stop, include the linear constraint that each stop is on exactly two trips. This means one arrival and one departure from the stop.

idxs = nchoosek(1:nStops,2);
dist = hypot(stopsLat(idxs(:,1)) - stopsLat(idxs(:,2)), ...
             stopsLon(idxs(:,1)) - stopsLon(idxs(:,2)));
lendist = length(dist);
Aeq = spones(1:length(idxs)); % Adds up the number of trips
beq = nStops;

Aeq = [Aeq;spalloc(nStops,length(idxs),nStops*(nStops-1))]; % allocate a sparse matrix
for ii = 1:nStops
    whichIdxs = (idxs == ii); % find the trips that include stop ii
    whichIdxs = sparse(sum(whichIdxs,2)); % include trips where ii is at either end
    Aeq(ii+1,:) = whichIdxs'; % include in the constraint matrix
end
beq = [beq; 2*ones(nStops,1)];

intcon = 1:lendist;
lb = zeros(lendist,1);
ub = ones(lendist,1);

opts = optimoptions('intlinprog','Display','off');
[xopt,costopt,exitflag,output] = intlinprog(dist,intcon,[],[],Aeq,beq,lb,ub,opts);

hold on
segments = find(xopt); % Get indices of lines on optimal path
lh = zeros(nStops,1); % Use to store handles to lines on plot
lh = updateSalesmanPlot(lh,xopt,idxs,stopsLon,stopsLat);
title('Solution with Subtours');

tours = detectSubtours(xopt,idxs);
numtours = length(tours); % number of subtours
fprintf('# of subtours: %d\n',numtours);

A = spalloc(0,lendist,0); % Allocate a sparse linear inequality constraint matrix
b = [];
while numtours > 1 % repeat until there is just one subtour
    % Add the subtour constraints
    b = [b;zeros(numtours,1)]; % allocate b
    A = [A;spalloc(numtours,lendist,nStops)]; % a guess at how many nonzeros to allocate
    for ii = 1:numtours
        rowIdx = size(A,1)+1; % Counter for indexing
        subTourIdx = tours{ii}; % Extract the current subtour
%         The next lines find all of the variables associated with the
%         particular subtour, then add an inequality constraint to prohibit
%         that subtour and all subtours that use those stops.
        variations = nchoosek(1:length(subTourIdx),2);
        for jj = 1:length(variations)
            whichVar = (sum(idxs==subTourIdx(variations(jj,1)),2)) & ...
                       (sum(idxs==subTourIdx(variations(jj,2)),2));
            A(rowIdx,whichVar) = 1;
        end
        b(rowIdx) = length(subTourIdx)-1; % One less trip than subtour stops
    end

    % Try to optimize again
    [xopt,costopt,exitflag,output] = intlinprog(dist,intcon,A,b,Aeq,beq,lb,ub,opts);

    % Visualize result
    lh = updateSalesmanPlot(lh,xopt,idxs,stopsLon,stopsLat);

    % How many subtours this time?
    tours = detectSubtours(xopt,idxs);
    numtours = length(tours); % number of subtours
    fprintf('# of subtours: %d\n',numtours);
end

title('Solusi mirza rofik'); 
hold off

