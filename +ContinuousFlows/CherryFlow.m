classdef CherryFlow < ContinuousFlows.AbstractODEFlow2D
%CHERRYFLOW Cherry Flow as described in Zaks (2001)
%
% dx = 1 + sin(y) - sin(x) sin(y)
% dy = b + cos(x) - K cos(x) cos(y)
%
%
% Standard domain: [0,2*pi] x [0, 2*pi] (as a 2-torus)
%
% b is the rotation number when the flow is Hamiltonian (K=1), and
% K the orientation of the flow (e.g., differentiating between source and sink)
%
% Parameter sets in constructor available:
% 'hamiltonian'
% 'sourcefractal'
% 'sinkfractal'
%

  properties
    % parameters for the flow
    K
    b
  end

  methods

    function obj = CherryFlow( dt, varargin )
    %CherryFlow Construct the cherry flow
    %
    % CherryFlow( dt, 'b', B, 'k', K )
    %
    % Construct Cherry flow with parameters
    % b (rotation number at K=1),
    % K (energy conservation parameter:
    %   -- K < 1 sink
    %   -- K = 1 center
    %   -- K > 1 source
    %
    % CherryFlow( dt, 'parset', PARAMETERSETNAME )
    % Params is:
    %
    % -- coefficient mu

        parser = inputParser;
        parser.addRequired('dt', @(x)isscalar(x) & (x > 0));
        parser.addOptional('parset',[], @ischar);
        parser.addParameter('b',(sqrt(5)-1)/2);
        parser.addParameter('K',1);

        parser.parse(dt,varargin{:});

      obj.dt = dt;
      obj.Domain = [0, 2*pi; 0,2*pi];

      % valid sets of parameters
      validsets.hamiltonian = {1, (sqrt(5)-1)/2};
      validsets.sourcefractal = {3/2,0.58766578};
      validsets.sinkfractal = {1/2,0.60398640};

      % try process the parameter set if possible
      try
          names = fieldnames(validsets);
          parset = validatestring( parser.Results.parset,...
              names );

          fprintf(1,'Cherry flow: %s set used\n', parset);
          myset = validsets.(parset);
          [obj.K, obj.b] = deal( myset{:} );

      catch e
          % if parameter set is not valid, use provided values
          if ~strcmpi(e.identifier,'MATLAB:unrecognizedStringChoice')
              rethrow(e)
          end
          disp('Valid set names: ')
          fieldnames(validsets)
          obj.b = parser.Results.b;
          obj.K = parser.Results.K;
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

        f(1,:) = 1 + sin(Y) - sin(X).*sin(Y);
        f(2,:) = obj.b + cos(X) - obj.K*cos(X).*cos(Y);

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


      J(1,1,:) =  -cos(X).*sin(Y);
      J(1,2,:) =   cos(Y) - sin(X).*cos(Y);

      J(2,1,:) =  -sin(X)+obj.K*sin(X).*cos(Y);
      J(2,2,:) =  obj.K*cos(X).*sin(Y);

      if nargin == 4
        Jj = jacobian@ContinuousFlows.AbstractODEFlow(obj, t, x, delta);
        DeltaJ = Jj - J;
      end


    end

  end

end
