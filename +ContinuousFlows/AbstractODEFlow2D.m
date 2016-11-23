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
    % DIVERGENCEPLOT(obj, t, R)
    %   As above, uses R points per axis of the obj.Domain (default: R =
    %   20).
    % DIVERGENCEPLOT(obj, t, xi, yi)
    %   As above, uses a tensor grid xi XX yi to plot.
    % [h] = DIVERGENCEPLOT(...)
    %   As above, returns graphics handle.
    % [X,Y,DIV] = DIVERGENCEPLOT(...)
    %   Returns spatial points and values of the divergence function.
    %   DIV is a matrix of size [rows(X), cols(X), numel(t)]
    %

    % compute grid based on input values
      if isempty(varargin)
        R = 100;
      elseif numel(varargin) == 1
        R = varargin{1};
      end

      if numel(varargin) < 2
        xi = linspace(obj.Domain(1,1), obj.Domain(1,2), R);
        yi = linspace(obj.Domain(2,1), obj.Domain(2,2), R);
      else
        assert( numel(varargin) == 2, 'We can use at most 4 arguments');
        xi = varargin{3};
        yi = varargin{4};
      end

      [X,Y] = meshgrid(xi, yi);

      x = [X(:),Y(:)].';

      Divs = nan( [size(X), numel(t)] );
      Div_i = nan( [size(X), 1] );

      for k = 1:numel(t)
        J = obj.jacobian(t(k),x);
        for xi = 1:size(x,2)
          Div_i(xi) = trace( J(:,:,xi) );
        end
        Divs(:,:,k) = reshape(Div_i,size(X));
      end

      if nargout > 1
        varargout = {X,Y,Divs};
      else
        for k = 1:numel(t)
          if k == 1
            [~,h] = contourf(X,Y,Divs(:,:,1));
          else
            h.Visible ='off';
            h.ZData = Divs(:,:,k);
            h.Visible = 'on';
          end
          title(sprintf('t = %.2f',t(k)));
          pause(1/15);
        end
        if nargout > 0
          varargout = h;
        end
        if nargout > 0
          varargout = h;
        end
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
    end

      [X,Y] = meshgrid(xi, yi);

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
    %
    % h = POLARVFPLOT(obj, t, R)
    %   As above, uses R points per axis of the obj.Domain (default: R =
    %   20).
    % h = POLARVFPLOT(obj, t, xi, yi)
    %   As above, uses a tensor grid xi XX yi to plot.
    % [h] = POLARVFPLOT(...)
    %   As above, returns graphics handle.
    % [X,Y,ANGLE,NORM] = POLARVFPLOT(...)
    %   Returns spatial points and values of the angle.
    %   OMEGA is a matrix of size [rows(X), cols(X), numel(t)]
    %


    % compute grid based on input values
      if isempty(varargin)
        R = 100;
      elseif numel(varargin) == 1
        R = varargin{1};
      end

      if numel(varargin) < 2
        xi = linspace(obj.Domain(1,1), obj.Domain(1,2), R);
        yi = linspace(obj.Domain(2,1), obj.Domain(2,2), R);
      else
        assert( numel(varargin) == 2, 'We can use at most 4 arguments');
        xi = varargin{3};
        yi = varargin{4};
      end

      [X,Y] = meshgrid(xi, yi);
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


  end % methods
end % classdef
