%% Using Continuous Flows
% ContinuousFlows is a MATLAB package of classes that simplify simulation of 
% dynamical systems, primarily flows derived from ODEs (but not exclusively).
% 
% This document is a demo of some of the basic features for new users.
% 
% % 
% you should add the folder *containing* the folder +ContinuousFlows into the 
% path

addpath('.');
%% 
% Dynamical systems are programmed as objects, on which you can perform operations 
% - creating visualizations/graphs, simulating, running diagnostics, etc. Since 
% many of you may be familiar with the Shadden--Marsden Double Gyre flow (2005) 
% let's use that as our example.
% 
% To read more about the model, you can run the help on it

help ContinuousFlows/DoubleGyre

%% 
% To read how to create the object, you should run help on the so-called "constructor 
% function", which means repeating the name of the model twice


help ContinuousFlows/DoubleGyre/DoubleGyre
%% 
% Some flows have a few parameters, other have many, including shorthands for 
% creating "sets" of parameters that are standard in literature. For example

dt = 0.01;
dynsys = ContinuousFlows.DoubleGyre(dt, 'standard');
%% 
% is the same as

dynsys = ContinuousFlows.DoubleGyre(dt, [0.1, 2*pi, 0.25]);

%% 
% The first argument is always "dt" - a time step at which the flow is simulated. 
% This can always be changed later by saying


dynsys.dt = 0.1; % for example

%% 
% If you want to see what parameters your model has, just call its name

dynsys
%% 
% Any of these parameters can be changed after the fact

dynsys.label = 'First Flow';
%% Plotting
% There are many plots that you can generate with a single command. To see the 
% range of functions available, call MATLAB's "methods" function

methods(dynsys)
%% 
% Then if you'd like to know more, call "help" on the specific function

help ContinuousFlows/DoubleGyre/quiverplot
%% 
% % 
% % 
% % 
% As an example, let's plot the velocity field using quiverplot. Notice that 
% help mentions "obj" as the first argument. This means you can call the function 
% in two formats:
%% 
% * quiverplot( dynsys, 0.1)
% * dynsys.quiverplot( 0.1 )
%% 
% % 
% Personally, I prefer the second, but it's really down to a preference.


figure(1)
dynsys.quiverplot(0.1); % at time = 0.1
colorbar;

% we can increase the resolution by changing the R parameter
% and create a video by requesting several timesteps at the same time
close(1)
dynsys.quiverplot(0:0.1:1, 40)
%% 
% Function I most commonly use is "flow": it simulates trajectories from (a 
% set of) initial conditions.

help dynsys.flow
%% 
% Let's use it to create a single trajectory

ic = [0.1; 0.2]; % column of coordinate values
T = 10; % endtime

[x,t] = dynsys.flow( ic, T ); % use t0 as the third parameter if your initial condition is not at t=0

%% 
% Trajectories are returned so that each *column* is a different timestep - 
% that's why we transpose each row before we plot (MATLAB uses rows as different 
% points in the plot)

hold on;
plot( x(1,:)', x(2,:)');
plot( x(1,1), x(2,1), 'ko'); % initial condition 
xlim(dynsys.Domain(1,:)); % set the plot limits based on the Domain parameter
ylim(dynsys.Domain(2,:));
hold off
%% 
% Often in our research we don't just simulate one or two trajectories, but 
% perhaps tens or hundreds. This package can do this efficiently - both initializing 
% a set of points on a grid or at random, and then simulating them.

figure(2);

%% 
% N points in a diamond shaped region


N = 10; % N points
diamond = [0.5,0.1; 0.9, 0.5; 0.5, 0.9; 0.1, 0.5];
source = polyshape( diamond(:,1), diamond(:,2) );
diamondIC = dynsys.samplePolygonInterior( N, diamond' );

%% 
% Plotting the points


diamond_shape= plot( polyshape( diamond(:,1), diamond(:,2) ),...
    'FaceColor','none', 'EdgeColor','b' ); hold on;
plot( diamondIC(1,:)',diamondIC(2,:)','db',...
    'MarkerFaceColor','b'); hold off;

%% 
% Simulate trajectories

Td = 1;
[xd,td] = dynsys.flow( diamondIC, Td );

%% 
% Now, let's see what the dimensions of xd are:

size(xd)
%% 
% Before, we simulated one i.c., so output was a matrix of 2 coord x K timesteps. 
% Now, the output is a 3D matrix, with third coordinate indexing different trajectories. 
% So, xd(:,:,5) corresponds to the trajectory of the 5th initial condition. We 
% can plot all of them together with a bit of MATLAB magic:


hold on;
plot( ...
    squeeze( xd(1,:,:) ), ... % extract all 1st coordinates into a flat matrix
    squeeze( xd(2,:,:) ), ... % extract all 2nd coordinates into a flat matrix
    'LineWidth',2);
hold off;
%% 
% There are other functions that you can use to generate initial conditions 
% - in the following list, they start with "sample":

methods(dynsys)

% for example, we can sample the domain at random (this works for systems
% that are in any number of dimensions):

domainIC = dynsys.sampleDomainGrid( 5, ... % 5 ^ dim initial conditions from a grid
    [0.1,1.9; ... % horizontal range
     0.1,0.9] );  % vertical range
hold on;
plot(domainIC(1,:)',domainIC(2,:)','sk','MarkerFaceColor','r');
[xx, tx] = dynsys.flow( domainIC, 2 );
plot( squeeze(xx(1,:,:)), squeeze(xx(2,:,:)), 'LineWidth',1);
hold off;
%% 
% There are other useful visualizations, particularly for 2D flows. They are 
% used very similarly to "quiverplot" visualization.

figure(3)

subplot(2,2,1);
dynsys.quiverplot(0); title('Velocity field');
subplot(2,2,2);
dynsys.streamplot(0); title('Stream function');
subplot(2,2,3);
dynsys.vorticityplot(0); title('Vorticity');
subplot(2,2,4);
dynsys.divergenceplot(0); title('Divergence');
%% 
% To create your own model, it's best to look at code of file  SimplestSystem.m 
% file and see how it's structured. 

simple = ContinuousFlows.SimplestSystem(0.1, 1 );

figure(4)
subplot(1,2,1);
simple.quiverplot(0); title('Velocity field');
subplot(1,2,2);
simple.vorticityplot(0); title('Vorticity');
%% 
% When you get familiar with the SimplestSystem.m, create a copy and edit that 
% file. Don't worry - you cannot break anything. In the worst case,you will re-download 
% ContinuousFlows from internet and have a fresh copy.
% 
% If you're working with 2D incompressible (divergence-free) system, then follow 
% the example of DoubleGyre.m.