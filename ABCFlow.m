%ABCFLOW
%  Unsteady ABC flow.
%
% dx = ( A + D t sin( pi t ) ) sin(z) + C cos(y)
% dy = B sin(x) + ( A + D t sin( pi t ) ) cos(z)
% dz = C sin(y) + B cos(x)

classdef ABCFlow < ContinuousFlows.ODEFlow

  properties
    % parameters for the flow
    A
    B
    C
    D
  end

  methods

    function obj = ABCFlow( dt, params )
    %ABCFLOW Construct an ABC flow object.
    % ABCFlow( dt, params )
    % Params can be:
    %
    % -- 1 x 4 vector of coefficients [A,B,C,D]
    % -- 'integrable' - integrable parameter set
    % -- 'steady' - chaotic autonomous parameter set
    % -- 'unsteady' chaotic non-autonomous parameter set

      obj.dt = dt;

      if ischar( params )
        switch params
          case 'integrable'
            params = [2, 1, 0, 0];
          case 'steady'
            params = [sqrt(3), sqrt(2), 1, 0];
          case 'unsteady'
            params = [sqrt(3), sqrt(2), 1, 1];
          otherwise
            error('Unknown parameter set');
        end
      end

      params = num2cell(params);
      [obj.A, obj.B, obj.C, obj.D] = deal(params{:});

      %% Set up integration parameters
      obj.integrator = @ode23t;
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

      tcoeff = obj.A + obj.D*t.*sin(pi*t);
      f(1,:) = tcoeff .* sin(x(3,:)) + obj.C .* cos(x(2,:));
      f(2,:) = tcoeff .* cos(x(3,:)) + obj.B .* sin(x(1,:));
      f(3,:) = obj.C*sin(x(2,:)) + obj.B*cos(x(1,:));

      % periodicity on the [0,1]^3 cube
      %      f = f/(2*pi);

    end

    function [ J ] = jacobian( obj, t, x )
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

      assert( numel(size(x)) == 2 );
      L = size(x,2);
      assert( numel(t) == 1 || numel(t) == L, ...
              ['Time is either a scalar or'...
               'has to match number of steps'] );

      tcoeff = obj.A + obj.D*t.*sin(pi*t);

      J(1,1,:) =  zeros(1,1,L);
      J(1,2,:) = -obj.C * sin( x(2,:) );
      J(1,3,:) =  tcoeff .* cos( x(3,:) );

      J(2,1,:) =  obj.B * cos( x(1,:) );
      J(2,2,:) =  zeros(1,1,L);
      J(2,3,:) = -tcoeff .* sin( x(3,:) );

      J(3,1,:) = -obj.B * sin( x(1,:) );
      J(3,2,:) =  obj.C * cos( x(2,:) );
      J(3,3,:) =  zeros(1,1,L);

    end

    function [varargout] = quiver( obj, t, N )
    %QUIVER Vector field of the flow.
    %
    % Produce vector field of the flow on the grid of N x N x N points at time
    % t.
    %
    % QUIVER(obj, t, N)
    %   Plots the vector field at time t on grid of N x N x N points.
    % h = QUIVER(obj, t, N)
    %   As above, and returns graphics handle of the quiver object.
    % [X,Y,Z,U,V,W] = QUIVER(obj, t, N)
    %   Returns spatial points and components of the vector field.
    %

      [X,Y,Z] = meshgrid(linspace(0,2*pi,N));
      f = obj.vf(t, [X(:),Y(:), Z(:)].');
      U = reshape(f(1,:), size(X));
      V = reshape(f(2,:), size(Y));
      W = reshape(f(3,:), size(Y));

      if nargout > 1
        varargout = {X,Y,Z,U,V,W};
      else
        h = quiver3(X,Y,Z,U,V,W);
        if nargout > 0
          varargout = h;
        end
      end

    end


  end

end
