classdef InterpolatedODEFlow2D < ContinuousFlows.AbstractODEFlow2D
%INTERPOLATEDODEFLOW2D Interpolated 2d velocity field, for example from PIV
%

  properties

    gridaxes % cell array of axis vectors
    Ux      % interpolated x-direction of the velocity field
    Uy      % interpolated y-direction of the velocity field
    isAutonomous % set to true if flow is autonomous (otherwise, last
                 % coordinate is time)
    period   % period of the time vector
    phase    % phase of the time vector, it is recorded as the value of the
             % first element in the time vector

  end

  methods

    function [varargout] = streamplot( obj, t, varargin)
    %STREAM Level sets of the stream function of the flow.
    %
    % STREAMPLOT(obj, t)
    %   Plots the streamlines of the interpolated flow starting from a
    %   random subset of gridpoints. If t has multiple elements, video is produced.
    % STREAMPLOT(obj,t, PARAMETER,VALUE )
    %   Additionally modify the behavior using named parameters. Default
    %   value in parenthesis
    %   'R' -- number of points in vertical/horizontal grid (100)
    %   'xi' -- manually set x grid; ignore if empty ([])
    %   'yi' -- manually set y grid; ignore if empty ([])
    %   'nlines' -- number of streamlines to compute (30)
    %   'stepsize' -- stepsize for streamline calcluation (see built-in
    %   STREAMLINE)
    %   'npoints' -- number of points in each streamline (see built-in STREAMLINE)
    %
    % [h] = STREAMPLOT(...)
    %   As above, returns graphics handle.
    % [X,Y,U,V] = STREAMPLOT(...)
    %   Returns spatial points and values of the velocity field.
    %   U,V are matrices of size [rows(X), cols(X), numel(t)]
    %
    % See also: STREAMLINE

      parser = inputParser;

      %% Process optional parameters
      parser.addParameter('R',100); % resolution of the grid
      parser.addParameter('xi', [] ); % xgrid that should be used
      parser.addParameter('yi', [] ); % ygrid that should be used
      parser.addParameter('nlines',30, @(n)(n>0));
      parser.addParameter('stepsize',0.1, @(x)(x>0 && x<1));
      parser.addParameter('npoints',200, @(n)(n>0));

      parser.parse(varargin{:});

      R = ceil(parser.Results.R);
      if isempty(parser.Results.xi)
        xi = linspace(obj.Domain(1,1), obj.Domain(1,2), R);
      else
        xi = parser.Results.xi;
      end
      if isempty(parser.Results.yi)
        yi = linspace(obj.Domain(2,1), obj.Domain(2,2), R);
      else
        yi = parser.Results.yi;
      end

      %% Compute 2D grid of velocities

      [X,Y] = meshgrid(xi, yi);
      x = [X(:),Y(:)].';
      U = nan( [size(X), numel(t)] );
      V = nan( [size(X), numel(t)] );
      for k = 1:numel(t)
        VF_i = obj.vf(t(k),x);
        U(:,:,k) = reshape(VF_i(1,:), size(X));
        V(:,:,k) = reshape(VF_i(2,:), size(X));
      end

      %% Use built-in streamline to construct streamlines

      if nargout > 1
        varargout = {X,Y,U, V};
      else
        sel = unique(randi(numel(X),ceil(parser.Results.nlines)));
        for k = 1:numel(t)
          cla;
          h = streamline(X,Y,...
                         U(:,:,k), V(:,:,k),...
                         X(sel), Y(sel), ...
                         [parser.Results.stepsize,...
                          ceil(parser.Results.npoints)]);
          title(sprintf('t = %.2f',t(k)));
          pause(1/15);
        end
        if nargout > 0
          varargout = {h};
        end
      end
    end


    function obj = InterpolatedODEFlow2D( dt, ...
                                          axesNodes, ...
                                          UxGrid, UyGrid, ...
                                          varargin )
      %INTERPOLATEDODEFLOW2D Construct interpolated velocity flow.
      %
      % obj = InterpolatedODEFlow2D( dt, ...
      %                              axesNodes, ...
      %                              UxGrid, UyGrid, ...
      %                              PROPERTY, VALUE, ... )
      % Arguments:
      % dt - timestep for the integrator
      % axesNodes - cell array of nodes along each axis (x, y, t)
      %           - if t is omitted, the flow is assumed autonomous
      % UxGrid    - x-velocity matrix array of size numel(x) X numel(y) [ X numel(t) ]
      % UyGrid    - y-velocity matrix array of size numel(x) X numel(y) [ X
      % numel(t) ]
      %
      % UxGrid and UyGrid follow the NDGRID standard (not MESHGRID!)
      %
      % Optional PROPERTY, VALUE pairs can be used to specify:
      % 'period'    - if specified, the input PIV data is taken to be
      %               periodic with the given period. This could be used
      %               if PIV samples only a single period of a periodic
      %               flow.
      %
      % 'interpolation' - interpolation method used by
      %                   griddedInterpolant. Please see griddedInterpolant
      %                   documentation for description of options (default:
      %                   'linear')
      %
      % 'extrapolation' - interpolation method used by
      %                   griddedInterpolant. Please see griddedInterpolant
      %                   documentation for description of options (default:
      %                   'nearest')
      %
      % 'test'          - plot the first slice of the quiver plot (with a
      %                   small grid offset) to visually compare input data
      %                   and interpolated data.
      %                   Passed value is the percentage of the grid step
      %                   which is used to offset data
      %                   (default: 0)
      %
      % The data should be ordered in such a manner that
      % quiver( axesNodes{1}, axesNodes{2}, UxGrid(:,:,1), UyGrid(:,:,2) )
      % produces a valid velocity field portrait.

      obj.dt = dt;
      validateattributes( UxGrid, {'numeric'}, {}, ...
                          'InterpolatedODEFlow2D', 'UxGrid', 3 );
      D = numel(size(UxGrid));

      validateattributes( axesNodes, {'cell'}, {'numel',D},...
                          'InterpolatedODEFlow2D', 'axesNodes', 2 );


      validateattributes( UyGrid, {'numeric'}, {'ndims',D},...
                          'InterpolatedODEFlow2D', 'UyGrid', 4 );

      assert( all(size(UxGrid) == size(UyGrid)),...
              'Velocity components have to be of the same size' );

      obj.isAutonomous = (D==2);

      parsearg = inputParser;
      parsearg.addParameter('period',NaN,@isnumeric);
      parsearg.addParameter('interpolation','linear',@ischar);
      parsearg.addParameter('extrapolation','nearest',@ischar);
      parsearg.addParameter('test',0,@(x)(x >= 0));

      parsearg.parse(varargin{:});
      params = parsearg.Results;

      obj.period = params.period;
      if ~isnan(obj.period)
        assert(~obj.isAutonomous, 'Periodization requires non-autnomous v.f.')
        obj.phase = axesNodes{3}(1); % phase is the first time step
      end

      % establish the domain
      obj.Domain = nan( 2, 2 );
      for d = 1:2
        obj.Domain(d,1) = min(axesNodes{d});
        obj.Domain(d,2) = max(axesNodes{d});
      end

      gridCols = size(UxGrid,2);
      gridRows = size(UxGrid,1);

      sizeX = numel(axesNodes{1});
      sizeY = numel(axesNodes{2});
      sizeT = numel(axesNodes{3});

      assert( sizeX == gridRows, 'X-coordinate should correspond to rows' );
      assert( sizeY == gridCols, 'Y-coordinate should correspond to columns' );
      assert( all(size(UxGrid) == size(UyGrid)),...
              'Dimensions of velocity matrices should match');

      if ~obj.isAutonomous
        assert( numel(axesNodes{3}) == size(UxGrid,3) );
      end

      % interpolate the velocity field
      obj.Ux = griddedInterpolant(axesNodes, ...
                                  UxGrid,...
                                  params.interpolation,...
                                  params.extrapolation);

      obj.Uy = griddedInterpolant(axesNodes, UyGrid,...
                                  params.interpolation,...
                                  params.extrapolation);

      %% Set up integration parameters
      obj.integrator = @ode23t;
      obj.intprops = odeset;
      obj.intprops = odeset(obj.intprops, 'Vectorized', 'on');
      %      obj.intprops = odeset(obj.intprops, 'Jacobian', @(t,x)obj.jacobian(t,x) );
      %      obj.intprops = odeset(obj.intprops, 'Stats','on' );

      if params.test > 0
        params.test = params.test/100;
        x = axesNodes{1};
        y = axesNodes{2};
        t = axesNodes{3};
        t0 = t(1);
        ts = t0+params.test*diff(t(1:2));


        [X,Y] = ndgrid(x,y);

        quiver(X,Y,UxGrid(:,:,1),UyGrid(:,:,1),'Color','r');
        hold all;

        h = obj.quiverplot(ts, 'grid', ...
                           {x+params.test*diff(x(1:2));...
                            y+params.test*diff(y(1:2));});
        h.Color = 'b';
        hold off;

        xlabel( sprintf('X - Range [%f,%f]', min(x), max(x) ) );
        ylabel( sprintf('Y - Range [%f,%f]', min(y), max(y) ) );
        title({sprintf('Input data (red) at t=%f',t0),...
               sprintf('Interpolated data (blue) at t = %f',ts)});

      end


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
        % periodize and take phase into account
        tt = rem( t - obj.phase, obj.period ) + obj.phase;

        % because of floating point precision,
        % tt may be different from ti even when rem has no effect
        % therefore check manually if the change is significant
        % and change the time input only when it is
        sel = abs(tt-t)/abs(obj.period) > 1e-8;
        t(sel) = tt(sel);
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

  end

end
