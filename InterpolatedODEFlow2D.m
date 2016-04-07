classdef InterpolatedODEFlow2D < ContinuousFlows.AbstractODEFlow2D
%INTERPOLATEDODEFLOW2D Interpolated 2d velocity field, for example from PIV
%

  properties

    gridaxes % cell array of axis vectors
    Ux      % interpolated x-direction of the velocity field
    Uy      % interpolated y-direction of the velocity field
    isAutonomous % set to true if flow is autonomous (otherwise, last
                 % coordinate is time)
    period   %

  end

  methods

    function obj = InterpolatedODEFlow2D( dt, ...
                                          axesNodes, ...
                                          UxGrid, UyGrid, ...
                                          period )
      %INTERPOLATEDODEFLOW2D Construct interpolated velocity flow.
      %
      % Arguments:
      % dt - timestep for the integrator
      % axesNodes - cell array of nodes along each axis (x, y, t)
      %           - if t is omitted, the flow is assumed autonomous
      % UxGrid    - x-velocity matrix array of size numel(x) X numel(y) [ X numel(t) ]
      % UyGrid    - y-velocity matrix array of size numel(x) X numel(y) [ X numel(t) ]
      % period    - (optional) desired period for the nonautonomous flow, e.g.,
      %             when the samples are within one interval of the velocity
      %             field period
      %
      %

      obj.dt = dt;

      validateattributes( axesNodes, {'cell'}, {} );
      D = numel(axesNodes);
      validateattributes( UxGrid, {'numeric'}, {'ndims',D} );
      validateattributes( UyGrid, {'numeric'}, {'ndims',D} );

      assert( all(size(UxGrid) == size(UyGrid)),...
              'Velocity components have to be of the same size' );

      obj.isAutonomous = (D==2);

      if exist('period','var') && ~isempty(period)
        obj.period = period;
      else
        obj.period = NaN;
      end

      if ~obj.isAutonomous
        obj.period = range(axesNodes{3});
      end

      % establish the domain
      obj.Domain = nan( 2, 2 );
      for d = 1:2
        obj.Domain(d,1) = min(axesNodes{d});
        obj.Domain(d,2) = max(axesNodes{d});
      end

      % interpolate the velocity field
      obj.Ux = griddedInterpolant(axesNodes, UxGrid);
      obj.Uy = griddedInterpolant(axesNodes, UyGrid);

      %% Set up integration parameters
      obj.integrator = @ode23t;
      obj.intprops = odeset;
      obj.intprops = odeset(obj.intprops, 'Vectorized', 'on');
      %      obj.intprops = odeset(obj.intprops, 'Jacobian', @(t,x)obj.jacobian(t,x) );
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

    % periodize time if needed
      if ~isnan(obj.period)
        t = mod( t, obj.period );
      end

      % construct the appropriate input vector for the velocity field
      X = x(1,:)';
      Y = x(2,:)';

      if isscalar(t)
        t = repmat(t, size(X));
      end

      % evaluate the velocity field and
      % assemble the output
      Uxv = obj.Ux(X,Y,t);
      Uyv = obj.Uy(X,Y,t);
      f = [Uxv(:), Uyv(:)]';
    end

    function [j] = jacobian( obj, t, x)
    % DISABLED

      j = [];

    end

  end

end
