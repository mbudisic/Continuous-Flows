classdef DMDdemoFlow < ContinuousFlows.AbstractHamiltonian2DFlow

  properties

    B % magnitude of the background flow
    P % magnitude of the pair of vortices
    G % magnitude of the growing vortex
    D % magnitude of the decaying vortex

    Ob % complex frequency of bg flow
    Op % complex frequency of pair of vortices
    Og % complex frequency of the growing vortex
    Od % complex frequency of the decaying vortex

  end

  methods

    function obj = DMDdemoFlow( dt, params )
    % parameter structure should have the following fields
    % B % magnitude of the background flow
    % P % magnitude of the pair of vortices
    % G % magnitude of the growing vortex
    % D % magnitude of the decaying vortex

    % Ob % complex frequency of bg flow
    % Op % complex frequency of pair of vortices
    % Og % complex frequency of the growing vortex
    % Od % complex frequency of the decaying vortex

      obj.B = 1/2;
      obj.P = 1/4;
      obj.G = 1/8;
      obj.D = 1/16;

      obj.Ob = 0;
      obj.Op = -(log(2)/20) + (2j*pi/2);
      obj.Og =  (log(2)/100) + (2j*pi/4);
      obj.Od = -log(2)/5 + (2j*pi/4);


      if nargin < 1
        help ContinuousFlows.HackbornWeldonRotOsc.HackbornWeldonRotOsc
      end

      if nargin == 2
        for fn = fieldnames(params)'
          name = fn{1};
          obj.(name) = params.(name);
        end
      end

      obj.Domain = [0,5; -1,1];
      obj.dt = dt;


      obj.integrator = @ode15s;
      obj.intprops = odeset;
      obj.intprops = odeset(obj.intprops, 'Vectorized', 'on');
      obj.intprops = odeset(obj.intprops, 'Jacobian', 'off');
    end

    function [out] = Psi( obj, t, x, o )
    % PSI Compute the stream function or its derivatives along a
    % trajectory given by (t, x)
    % [ out ] = Psi( obj, t, x, order )
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


      out = obj.P*obj.Gauss( t, x, o, [1;-0.5], eye(2)*0.1, obj.Op ) + ...
            obj.P*obj.Gauss( t, x, o, [1; 0.5], eye(2)*0.1, conj(obj.Op) ) + ...
            obj.G*obj.Gauss( t, x, o, [2; 0], eye(2)*0.1, obj.Og ) + ...
            obj.D*obj.Gauss( t, x, o, [4; 0], eye(2)*0.1, obj.Od );

      switch o
        case 0
          out = out + obj.B*x(2,:);
        case 1
          out = out + obj.B*[0;1];
        otherwise
          error('Derivatives 0,1,2 are the only implemented');
      end

    end

    function [out] = Gauss( obj, t, x, o, mu, sigma, omega )


      D = size(x,1);
      assert( all( size(mu) == [D,1] , 'all' ));
      assert( all( size(sigma) == [D,D] ,'all' ));

      mu = reshape(mu, [],1);
      o0 = @(x)mvnpdf( x.', mu, sigma ).';

      % time evolution
      tcoeff = exp(omega*reshape(t,1,[]));
      if imag(omega) >= 0,
        tcoeff = real(tcoeff);
      else
        tcoeff = imag(tcoeff);
      end
      isigma = pinv(sigma);

      p = @(x)sqrt( exp(diag( -(x-mu).' * isigma * (x-mu) ) )...
                    /det(2*pi*sigma) ).';

      switch(o)
        case 0
          V = p(x);
        case 1
          V = -p(x).*( isigma *(x-mu) );
        otherwise
          error('Derivatives 0,1,2 are the only implemented');
      end

      out = V.*tcoeff;
      assert( size(out,2) == size(x,2) );

    end

  end

end