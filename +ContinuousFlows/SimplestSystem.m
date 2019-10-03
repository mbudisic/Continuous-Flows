classdef SimplestSystem < ContinuousFlows.AbstractODEFlow2D
%SIMPLESTSYSTEM A simple extension to demonstrate a 2D flow.
%
% dx = alpha*y
% dy = - alpha*x
%

  properties
    % parameters for the flow
    alpha
  end

  methods

    function obj = SimplestSystem( dt, alpha )
    %SimplestSystem Construct the system.
    % SimplestSystem( dt, params )
    % Params is:
    %
    % -- coefficient alpha

    % you don't need to modify this
      obj.dt = dt; % time step
      obj.Domain = [-6,6; -6,6]; % typical range for the state space
    
    % set this to parameters you need
      obj.alpha = alpha;

      %% Set up integration parameters - no need to change this
      obj.integrator = @ode45;

    end

    function [ f ] = vf( obj, t, x )
    % Compute vector field along
    % a single trajectory given by (t, x)
    %
    % this is the function you want to change to set up equations for your
    % system


      f(1,:) = obj.alpha * x(2,:);
      f(2,:) = -obj.alpha * x(1,:);

    end


  end

end
