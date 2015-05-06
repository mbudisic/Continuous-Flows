%DOUBLEGYRE
%  Planar double-gyre system, as used by Shadden (2005)
%

classdef DoubleGyre < ContinuousFlows.ODEFlow

  properties
    A
    omega
    epsilon
  end

  methods

    function obj = DoubleGyre( dt, params )
    %DOUBLEGYRE Construct a Double Gyre object.
    % DoubleGyre( dt, params )
    % Params can be
    %
    % -- 1 x 3 vector of coefficients [A,omega, epsilon]
    % -- 'standard' - parameter set [0.1, 2*pi, 0.25]

      obj.dt = dt;
      if ischar( params )
        switch params
          case 'standard'
            params = [0.1, 2*pi, 0.25]; % frequency == 1
          otherwise
            error('Unknown parameter set');
        end
      end

      params = num2cell(params);
      [obj.A, obj.omega, obj.epsilon] = deal(params{:});

      %% Set up integration parameters
      obj.integrator = @ode23t;
      obj.intprops = odeset;
      obj.intprops = odeset(obj.intprops, 'Vectorized', 'on');
      obj.intprops = odeset(obj.intprops, 'Jacobian', @(t,x)obj.jacobian(t,x) );
      %      obj.intprops = odeset(obj.intprops, 'Stats','on' );

    end

    function [ f ] = vf( obj, t, x )
    % VF Compute vector field along a trajectory
    % a single trajectory given by (t, x)
    % [ f ] = vf( obj, t, x )
    %
    % t   - row-vector of times
    % x   - trajectory
    %     - columns correspond to time steps
    %     - rows correspond to states
    %
    % Returns:
    % f   - evaluation of the vector field
    %     - each f(:,i) is a dim x 1 vector field evaluation
    %     - of the vector field at [ t(i), x(i,:) ] point

    %%
    % Time-varying coefficients
      a = obj.epsilon .* sin( obj.omega * t );
      b = 1 - 2 * a;

      F1 = a.*x(1,:).^2 + b.*x(1,:);
      F2 = 2*a.*x(1,:) + b;

      %%
      % Velocity field
      f(1,:) = -pi*obj.A*sin(pi*F1) .* cos(pi*x(2,:));
      f(2,:) =  pi*obj.A*cos(pi*F1) .* sin(pi*x(2,:)).*F2;
    end

    function [ J ] = jacobian( obj, t, x )
    % JACOBIAN Compute Jacobian of the vector field along
    % a single trajectory given by (t, x)
    % [ J ] = jacobian( obj, t, x )
    %
    % t   - row-vector of times
    % x   - trajectory
    %     - columns correspond to time steps
    %     - rows correspond to states
    % Returns:
    % J   - Jacobians
    %     - each J(:,:,i) is a dim x dim Jacobian matrix
    %     - of the vector field at [ t(i), x(i,:) ] point

      assert( ismatrix(x) );
      L = size(x,2);
      assert( numel(t) == 1 || numel(t) == L, ...
              ['Time is either a scalar or'...
               'has to match number of steps'] );
    %%
    % Time-varying coefficients
      a = obj.epsilon .* sin( obj.omega * t );
      b = 1 - 2 * a;

      F1 = a.*x(1,:).^2 + b.*x(1,:);
      F2 = 2*a.*x(1,:) + b;

      J(1,1,:) = (-pi^2*obj.A)*cos(pi*F1).*cos(pi*x(2,:) ).*F2;
      J(1,2,:) = (pi^2*obj.A)*sin(pi*F1).*sin(pi*x(2,:));
      J(2,1,:) = (-pi^2*obj.A)*sin(pi*F1).*sin(pi*x(2,:)).*F2.^2 ...
          + (2*pi*obj.A*a).*cos(pi*F1).*sin(pi*x(2,:));
      J(2,2,:) = (pi^2*obj.A).*cos(pi*F1).*cos(pi*x(2,:)).*F2;

    end

    function [varargout] = quiver( obj, t, N )
    %QUIVER Vector field of the flow.
    %
    % Produce vector field of the flow in 2*N x N points
    % at time t.
    %
    % QUIVER(obj, t, N)
    %   Plots the vector field at time t on grid of N x N points.
    % h = QUIVER(obj, t, N)
    %   As above, and returns graphics handle of the quiver object.
    % [X,Y,U,V] = QUIVER(obj, t, N)
    %   Returns spatial points and components of the vector field.


      [X,Y] = meshgrid(linspace(0,2,N), ...
                       linspace(0,1,N) );
      f = obj.vf(t, [X(:),Y(:)].');
      U = reshape(f(1,:), size(X));
      V = reshape(f(2,:), size(Y));

      if nargout > 1
        varargout = {X,Y,U,V};
      else
        h = quiver(X,Y,U,V);
        if nargout > 0
          varargout = h;
        end
      end

    end
  end

end
