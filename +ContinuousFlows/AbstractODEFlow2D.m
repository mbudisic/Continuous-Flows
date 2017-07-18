classdef (Abstract) AbstractODEFlow2D < ContinuousFlows.AbstractODEFlow
%VELOCITYFIELD2D Abstract methods for 2D plots of velocity field.

  methods

    function [varargout] = divergenceplot( obj, t, varargin)
    %DIVERGENCEPLOT Level sets of the divergence function of the flow.
    %
    % DIVERGENCEPLOT(obj, t)
    %   Plots the divergence at time t on the default grid on
    %   obj.Domain.
    %   If t has multiple elements, video is produced.
    % DIVERGENCEPLOT(...,'R', R)
    %   As above, uses R points per axis of the obj.Domain (default: R =
    %   20).
    % DIVERGENCEPLOT(...,'grid', {xi,yi} )
    %   As above, uses a tensor grid xi XX yi to plot. xi and yi are
    %   1D vectors.
    % [h] = DIVERGENCEPLOT(...,'normalized',true)
    %   Normalizes divergence by magnitude of velocity.
    % [X,Y,DIV] = DIVERGENCEPLOT(...)
    %   Returns spatial points and values of the divergence function.
    %   DIV is a matrix of size [rows(X), cols(X), numel(t)]
    %

      parser = inputParser;
      parser.addRequired('t');
      parser.addParameter('R',100, @(x)x>0);
      parser.addParameter('grid',{}, @iscell);
      parser.addParameter('normalized',false,@islogical);

      parser.parse(t, varargin{:});

      params = parser.Results;
      % compute grid based on input values

      if isempty(params.grid)
        xi = linspace(obj.Domain(1,1), obj.Domain(1,2), params.R);
        yi = linspace(obj.Domain(2,1), obj.Domain(2,2), params.R);
      else
        xi = params.grid{1};
        yi = params.grid{2};
        validateattributes( xi, {'numeric'},{'vector','real'});
        validateattributes( yi, {'numeric'},{'vector','real'});
      end

      [X,Y] = ndgrid(xi, yi);

      x = [X(:),Y(:)].';

      Divs = nan( [size(X), numel(t)] ); % divergences
      Div_i = nan( [size(X), 1] ); % unsorted divergences

      for k = 1:numel(t)
        J = obj.jacobian(t(k),x);
        for ii = 1:size(x,2)
          Div_i(ii) = trace( J(:,:,ii) );
        end
        Divs(:,:,k) = reshape(Div_i,size(X));

        % normalize by velocity magnitude if requested
        if (params.normalized)
          f = obj.vf(t(k),x);
          f_norm = reshape( hypot( f(1,:), f(2,:) ), size(X) );
          Divs(:,:,k) = Divs(:,:,k) ./ f_norm;
        end

        %% plotting
        if nargout <= 1
          if k == 1
            h = image('XData',xi,...
                      'YData',yi, ...
                      'CData',Divs(:,:,1),...
                      'CDataMapping','scaled');
            xlim([min(xi),max(xi)]);
            ylim([min(yi),max(yi)]);

            hb = colorbar;

            if (params.normalized)
              title(hb,{'Divergence/','Magnitude'})
            else
              title(hb,{'Divergence'})
            end

            divslice = abs(Divs(:,:,1));
            md = median(divslice(:));
            caxis([-md,md]);


          else
            h.Visible ='off';
            h.CData = Divs(:,:,k);
            h.Visible = 'on';

          end
          title(sprintf('t = %.2f',t(k)));
          pause(1/15);
        end

      end

      if nargout == 1
        varargout = h;
      elseif nargout > 1
        varargout = {X,Y,Divs};
      end
    end

    function [varargout] = quiverplot( obj, t, varargin )
    %QUIVERPLOT velocity field of the flow.
    %
    % Produce velocity field of the flow at time t on the tensor product grid
    % xi XX yi
    %
    % QUIVERPLOT(obj, t)
    %   Plots the velocity field at time t on the default grid on
    %   obj.Domain.
    %   If t has multiple elements, video is produced.
    % QUIVERPLOT(obj, t, R)
    %   As above, uses R points per axis of the obj.Domain (default: R =
    %   20).
    % QUIVERPLOT(..., 'grid', {xi, yi})
    %   As above, uses a tensor grid xi XX yi to plot (Ignores R parameter)
    % QUIVERPLOT(..., 'directionOnly', true)
    %   As above, but scale all vectors to unit length.
    %
    % [h] = QUIVERPLOT(...)
    %   As above, returns graphics handle.
    % [X,Y,U,V] = QUIVERPLOT(...)
    %   Returns spatial points and components of the velocity field.
    %   U and V are matrices of size [rows(X), cols(X), numel(t)]

      parser = inputParser;
      parser.addRequired('t');
      parser.addOptional('R',20, @isscalar);
      parser.addParameter('grid',{},@iscell);
      parser.addParameter('directionOnly',false,@islogical);

      parser.parse(t, varargin{:});
      params = parser.Results;

      R = params.R;

      if numel(params.grid) ~= 2
        xi = linspace(obj.Domain(1,1), obj.Domain(1,2), params.R);
        yi = linspace(obj.Domain(2,1), obj.Domain(2,2), R);
      else
        xi = params.grid{1};
        yi = params.grid{2};
        validateattributes( xi, {'numeric'},{'vector','real'});
        validateattributes( yi, {'numeric'},{'vector','real'});
      end

      [X,Y] = ndgrid(xi, yi); % use NDGrid format

      U = nan( [size(X), numel(t)] );
      V = nan( [size(X), numel(t)] );

      for k = 1:numel(t)
        f = obj.vf(t(k), [X(:),Y(:)].');
        U(:,:,k) = reshape(f(1,:), size(X));
        V(:,:,k) = reshape(f(2,:), size(Y));
      end

      if nargout > 1
        varargout = {X,Y,U,V};
      else
        for k = 1:numel(t)
          if k == 1
            h = quiver(X,Y,U(:,:,1),V(:,:,1));
            axis manual;
            %            h.AutoScale = 'off';
          else
            h.Visible = 'off';
            h.UData = U(:,:,k);
            h.VData = V(:,:,k);
            h.Visible = 'on';
          end
          if params.directionOnly
            M = hypot(h.UData, h.VData);
            h.UData = 0.8 * h.UData ./ M;
            h.VData = 0.8 * h.VData ./ M;
            %            h.AutoScale = 'off';
          end
          title(sprintf('t = %.2f',t(k)));
          pause(1/15);
        end
        if nargout > 0
          varargout = {h};
        end % if
      end % if
    end % quiverplot

    function [varargout] = polarvfplot( obj, t, varargin)
    %POLARVFPLOT Level sets of the angle of the velocity flow.
    %
    % POLARVFPLOT(obj, t)
    %   Plots the scalar (z-component) angle field at time t on the
    %   default grid on obj.Domain.
    %   If t has multiple elements, video is produced.
    % POLARVFPLOT(...,'R', R)
    %   As above, uses R points per axis of the obj.Domain (default: R =
    %   20).
    % POLARVFPLOT(...,'grid', {xi,yi} )
    %   As above, uses a tensor grid xi XX yi to plot. xi and yi are
    %   1D vectors.
    % [h] = POLARVFPLOT(...,'logarithmic',true)
    %   Plots color scale based on the log-magnitude of the velocity.
    % [X,Y,Omega,Norm] = POLARVFPLOT(...)
    %   Returns spatial points and velocity angle and magnitude (or its logarithm).
    %   Omega, Norm are matrices of size [rows(X), cols(X), numel(t)]
    %

      parser = inputParser;
      parser.addRequired('t');
      parser.addParameter('R',20, @(x)x>0);
      parser.addParameter('grid',{}, @iscell);
      parser.addParameter('logarithmic',false,@islogical);

      parser.parse(t, varargin{:});
      params = parser.Results;

      t = params.t;

      % compute grid based on input values

      if isempty(params.grid)
        xi = linspace(obj.Domain(1,1), obj.Domain(1,2), params.R);
        yi = linspace(obj.Domain(2,1), obj.Domain(2,2), params.R);
      else
        xi = params.grid{1};
        yi = params.grid{2};
        validateattributes( xi, {'numeric'},{'vector','real'});
        validateattributes( yi, {'numeric'},{'vector','real'});
      end

      [X,Y] = ndgrid(xi, yi);

      x = [X(:),Y(:)].';

      Omega = nan( [size(X), numel(t)] );
      Norm = nan( [size(X), numel(t)] );

      for k = 1:numel(t)
        v = obj.vf(t(k),x);
        w = complex(v(1,:),v(2,:));
        Omega_i = angle( w );
        Norm_i = abs(w);
        Omega(:,:,k) = reshape(Omega_i,size(X));
        Norm(:,:,k) = reshape(Norm_i,size(X));
      end

      if parser.Results.logarithmic
        Norm = log10(Norm);
      end

      if nargout > 1
        varargout = {X,Y,Omega,Norm};
      else
        for k = 1:numel(t)
          if k == 1
            [~,h.lines] = contour(X,Y,Norm(:,:,1),'k');
            hold on;
            h.color = pcolor(X,Y,Omega(:,:,1));
            h.color.ZData(:) = -0.1;
            shading flat;
            colormap(hsv(1024)); caxis([-pi,pi]);
            colorbar;
          else
            h.color.Visible ='off';
            h.color.CData = Omega(:,:,k);
            h.lines.ZData = squeeze(Norm(:,:,k));
            h.color.Visible = 'on';
          end
          title(sprintf('t = %.2f',t(k)));
          pause(1/15);
        end
        if nargout > 0
          varargout = {h};
        end
      end
    end % polarvfplot

    function Points = samplePolygonBoundary( obj, N, polygon )
    %SAMPLEPOLYGONBOUNDARY Return uniform sample of points along polygon
    %boundary of a domain.
    %
    % Points = obj.samplePolygonBoundary( N, polygon )
    % Returns a 2 x N matrix of points sampled along the boundary of the
    % polygon (2 x M matrix), uniformly gridded as possible.
    %
    %
    % See also: sampleDomainGrid

      validateattributes( polygon, ...
                          {'numeric'}, ...
                          {'nrows',2, 'nonnan','finite'},'',...
                          'polygon',3);

      S = size(polygon,2); % number of sides

      validateattributes( N, {'numeric'}, {'>=', size(polygon, 2)},'',...
                          'N',2 );

      % compute lengths of sides
      P = [polygon, polygon(:,1)];
      D = sum( (P(:,2:end) - P(:,1:end-1)).^2, 1 ).^(1/2);

      % compute number of points along each side
      ni = fix( N * D / sum(D) );

      % if points need to be subtracted, subtract them from the
      % smallest-gap side
      while sum(ni) > N
        [~,argmin] = min( D ./ ni );
        ni(argmin) = ni(argmin) - 1;
      end

      % if points need to be added, add them to the largest-gap side
      while sum(ni) < N
        [~,argmax] = max( D ./ ni );
        ni(argmax) = ni(argmax) + 1;
      end

      % linearly sample each side and add to the output
      Points = [];
      for s = 1:S
        SidePoints = [ linspace(P(1,s), P(1,s+1), ni(s)+1); ...
                       linspace(P(2,s), P(2,s+1), ni(s)+1) ];
        Points = [Points, SidePoints(:,1:end-1)];
      end

    end

    function Points = samplePolygonInterior( obj, N, polygon )
    %SAMPLEPOLYGONINTERIOR Return a sample of initial conditions inside the polygon.
    %
    % Points = obj.samplePolygonBoundary( N, polygon, samplefun )
    % Returns a 2 x N matrix of points sampled on the inside of the
    % polygon, using sample-and-reject technique.
    %
    % See also: samplePolygonaBoundary, sampleDomainRandom

      validateattributes( polygon, ...
                          {'numeric'}, ...
                          {'nrows',2, 'nonnan','finite'},'',...
                          'polygon',3);

      S = size(polygon,2); % number of sides

      domain = [min(polygon(1,:)), max(polygon(1,:)); ...
                min(polygon(2,:)), max(polygon(2,:)) ];

      PolyArea = polyarea( polygon(1,:), polygon(2,:) );
      DomainArea = range(domain(1,:))*range(domain(2,:));
      OversamplingRatio = DomainArea/PolyArea;

      Points = [];

      while size(Points,2) < N

        % sample the surrounding square domain
        NewSample = obj.sampleDomainRandom( fix(OversamplingRatio*N),...
                                            domain );

        % choose the points inside the polygon
        sel = inpolygon( NewSample(1,:), NewSample(2,:),...
                         polygon(1,:), polygon(2,:) );

        NewSample = NewSample(:, sel);

        % add the required number of points
        Deficit = N - size(Points,2);
        Points = [Points, NewSample(:,1:min(Deficit, size(NewSample,2)))];
      end
    end


  end % methods
end % classdef
