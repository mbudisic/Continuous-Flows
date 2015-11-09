classdef FourGyre < ContinuousFlows.Hamiltonian2DFlow
%FOURGYRE Planar four-gyre system, as used by Mezic (2010)
%
% Standard domain: [0,1] x [0,1]

  properties
    epsilon

  end

  methods

    function obj = FourGyre( dt, params )
    %DOUBLEGYRE Construct a Double Gyre object.
    % FourGyre( dt, params )
    % Params can be
    %
    % -- 1 x 3 vector of coefficients [A,omega, epsilon]
    % -- 'standard' - parameter set [0.1, 2*pi, 0.25]
    % -- 'steady'   - parameter set [0.1, 2*pi, 0.00]


      if nargin < 2
        help ContinuousFlows.FourGyre.FourGyre
      end

      obj.dt = dt;
      obj.Domain = [0,1; 0,1];

      if ischar( params )
        switch params
          otherwise
            error('Unknown parameter set');
        end
      end

      params = num2cell(params);
      [obj.epsilon] = deal(params{:});

      %% Set up integration parameters
      obj.integrator = @ode23t;
      obj.intprops = odeset;
      obj.intprops = odeset(obj.intprops, 'Vectorized', 'on');
      obj.intprops = odeset(obj.intprops, 'Jacobian', @(t,x)obj.jacobian(t,x) );
      %      obj.intprops = odeset(obj.intprops, 'Stats','on' );

    end

    function [out] = Psi( obj, t, x, o )
    % PSI Compute the stream function or its derivatives along a
    % trajectory given by (t, x)
    % [ f ] = Psi( obj, t, x, order )
    %
    % t   - row-vector of times
    % x   - trajectory
    %     - columns correspond to time steps
    %     - rows correspond to states
    % o   - order of calculation
    %
    % Returns:
    % out   - evaluation of the stream function or its derivatives
    %       - if o == 0, out is 1 x Nx row vector
    %       - if o == 1, out is 2 x Nx row vector; rows are x and y derivatives
    %                    respectively
    %       - if o == 2, out is 3 x Nx row vector; rows are xx, xy, yy
    %                    derivatives respectively

    %      Nx = size(x,2);

      X = x(1,:);
      Y = x(2,:);

      COSX = cos(2*pi*X);
      SINX = sin(2*pi*X);
      COSY = cos(2*pi*Y);
      SINY = sin(2*pi*Y);

      COSSIN = COSX.*SINY;
      SINSIN = SINX.*SINY;
      COSCOS = COSX.*COSY;
      SINCOS = SINX.*COSY;

      TC = obj.epsilon*cos(2*pi*t);

      if o == 0
        out = (SINSIN + TC.*COSCOS)/(2*pi);
      elseif o == 1
        dXPsi = COSSIN - TC.*SINCOS;
        dYPsi = SINCOS - TC.*COSSIN;
        out = [dXPsi; dYPsi];
      elseif o == 2
        dXXPsi = -SINSIN - TC.*COSCOS;
        dXYPsi =  COSCOS + TC.*SINSIN;
        dYYPsi = -SINSIN - TC.*COSCOS;
        out = 2*pi*[dXXPsi; dXYPsi; dYYPsi];
      else
        error('Higher orders not implemented');
      end
    end
  end

end
