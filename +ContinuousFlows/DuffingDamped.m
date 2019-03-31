classdef DuffingDamped < ContinuousFlows.AbstractODEFlow2D
%DUFFINGDAMPED Damped Duffing Oscillator
%
% dx = y
% dy = -delta * y - beta*x - alpha*x^3 + A sin(omega t)
%
%
% Standard domain: [-2,2] x [-1, 1]
%
%
% Parameter sets in constructor available:
% 'simple'
% 'mixing'
%

  properties
    % parameters for the flow
    alpha
    beta
    delta
    A
    omega
  end

  methods

    function obj = DuffingDamped( dt, varargin )
    %DuffingDamped Construct the duffing flow

        parser = inputParser;
        parser.addRequired('dt', @(x)isscalar(x) & (x > 0));
        parser.addOptional('parset',[], @ischar);
        parser.addParameter('alpha', 0);
        parser.addParameter('beta', 0);
        parser.addParameter('delta', 0);
        parser.addParameter('A', 0);
        parser.addParameter('omega', 0);

        parser.parse(dt,varargin{:});

      obj.dt = dt;
      obj.Domain = [-2,2;-1,1];

      % valid sets of parameters
      validsets.simple = {1,-1,0,0,0};
      validsets.mixing = {1,-1,0.5,0.42,1};

      % try process the parameter set if possible
      try
          names = fieldnames(validsets);
          parset = validatestring( parser.Results.parset,...
              names );

          fprintf(1,'Parameter set %s used\n', parset);
          myset = validsets.(parset);
          [obj.alpha, obj.beta, obj.delta, obj.A, obj.omega] = deal( myset{:} );

      catch e
          % if parameter set is not valid, use provided values
          if ~strcmpi(e.identifier,'MATLAB:unrecognizedStringChoice')
              rethrow(e)
          end
          disp('Valid set names: ')
          fieldnames(validsets)
          obj.alpha = parser.Results.alpha;
          obj.beta = parser.Results.beta;
          obj.delta = parser.Results.delta;
          obj.A = parser.Results.A;
          obj.omega = parser.Results.omega;
      end


      %% Set up integration parameters
      obj.integrator = @ode113;
      obj.intprops = odeset;
      obj.intprops = odeset(obj.intprops, 'Vectorized', 'on');
      obj.intprops = odeset(obj.intprops, 'Jacobian', @(t,x)obj.jacobian(t,x) );
      obj.intprops = odeset(obj.intprops, 'AbsTol',1e-12);
      obj.intprops = odeset(obj.intprops, 'RelTol',1e-6);
      %      obj.intprops = odeset(obj.intprops, 'Stats','on');

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

        X = x(1,:);
        Y = x(2,:);

        f(1,:) = Y;
        f(2,:) = -obj.delta*Y - obj.beta*X -obj.alpha*X.^3 + obj.A * ...
                 sin(obj.omega * t);

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

      X = x(1,:);
      Y = x(2,:);


      J(1,1,:) =  0;
      J(1,2,:) =  1;

      J(2,1,:) =  -obj.beta - 2*obj.alpha*X.^2;
      J(2,2,:) =  -obj.delta;

      if nargin == 4
        Jj = jacobian@ContinuousFlows.AbstractODEFlow(obj, t, x, delta);
        DeltaJ = Jj - J;
      end


    end

  end

end
