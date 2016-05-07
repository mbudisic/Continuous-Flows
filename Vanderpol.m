classdef Vanderpol < ContinuousFlows.AbstractODEFlow2D
%VANDERPOL Vanderpol oscillator.
%
% dx = y
% dy = mu ( 1 - x^2 )y - x
%

  properties
    % parameters for the flow
    mu
  end

  methods

    function obj = Vanderpol( dt, params )
    %Vanderpol Construct a Van der Pol system.
    % Vanderpol( dt, params )
    % Params is:
    %
    % -- coefficient mu

      obj.dt = dt;
      obj.Domain = [-6,6; -6,6];

      obj.mu = params;

      %% Set up integration parameters
      obj.integrator = @ode45;
      obj.intprops = odeset;
      obj.intprops = odeset(obj.intprops, 'Vectorized', 'on');
      obj.intprops = odeset(obj.intprops, 'Jacobian', @(t,x)obj.jacobian(t,x) );
      %      obj.intprops = odeset(obj.intprops, 'Stats','on' );

    end

    function [ f ] = vf( obj, t, x )
    % Compute vector field along
    % a single trajectory given by (t, x)
    %
    % t   - row-vector of times OR a scalar
    % x   - trajectory
    %     - rows correspond to time steps
    %     - columns correspond to states
    %
    % Returns:
    % f   - evaluation of the vector field
    %     - each f(:,i) is a dim x 1 vector field evaluation
    %     - of the vector field at [ t(i), x(i,:) ] point

    % periodicity on the [0,1]^3 cube
    %      x = 2*pi*x;

      f(1,:) = x(2,:);
      f(2,:) = obj.mu*(1 - x(1,:).^2).*x(2,:) - x(1,:);

    end

    function [ J, DeltaJ ] = jacobian( obj, t, x, delta )
    %JACOBIAN Compute Jacobian of the vector field along
    % a single trajectory given by (t, x)
    %
    % [ J ] = obj.jacobian( t, x )
    %
    % t   - row-vector of times OR a scalar
    % x   - trajectory
    %     - columns correspond to time steps
    %     - rows correspond to states
    % Returns:
    % J   - Jacobians
    %     - each J(:,:,i) is a dim x dim Jacobian matrix
    %     - of the vector field at [ t(i), x(i,:) ] point
    %
    % [ J, DeltaJ ] = obj.jacobian( ..., delta )
    %
    % Additionally return the difference between the numerically computed
    % and analytically computed Jacobians.
    %

      assert( numel(size(x)) == 2 );
      L = size(x,2);
      assert( numel(t) == 1 || numel(t) == L, ...
              ['Time is either a scalar or'...
               'has to match number of steps'] );

      J(1,1,:) =  zeros(1,1,L);
      J(1,2,:) =  ones(1,1,L);

      J(2,1,:) =  -2*obj.mu*x(1,:).*x(2,:) - 1;
      J(2,2,:) =  obj.mu.*(1-x(1,:).^2);

      if nargin == 4
        Jj = jacobian@ContinuousFlows.AbstractODEFlow(obj, t, x, delta);
        DeltaJ = Jj - J;
      end


    end

  end

end
