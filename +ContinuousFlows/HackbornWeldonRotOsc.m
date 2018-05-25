classdef HackbornWeldonRotOsc < ContinuousFlows.AbstractHamiltonian2DFlow
%HACKBORNWELDONROTOSC Hackborn-Weldon Rotor-Oscillator flow -- Hackborn et
%al., JFM, (1997), Weldon (2008) JFM
%
% The flow evolves in the channel [-1,1] x [-inf, inf] although
% practically [-1,1] x [-4, 4] is enough for common parameters. In this
% implementation, x and y variables are placed into state vector as
% [ y ]
% [ x ]
% in order to make the channel flow horizontally.
%
% Stream function $\Psi(x,y,t)$ is given by three components:
% $$ \Psi(x,y,t) = \Phi(x,y) + \Gamma(x,y) + \epsilon * \lambda * \Lambda(t,x) $$
%
% Free rotlet:
% $$ \Phi(x,y) = (1/2) \log \frac{\cosh(\pi y/2) - \cos[\pi(x-c)/2]}{\cosh(\pi
% y/2) + \cos[\pi(x+c)/2]} $$
%
% No-slip boundary correction to rotlet
% $$ \Gamma(x,y) = \int_0^\infty \cos(k y) G(x,k) $$
% where
% $$ G(x,k) = ... $$
%
% This integral is evaluated using Legendre-Gauss weights written by Greg von
% Winckel's (as lgwt)
%
% Wall-induced periodic shear
% $$\Lambda(t,x) = (A*x + B*x^2/2) \cos(\lambda t)
% sets up a background linear shear profile, that oscillates in magnitude.
%
% In Hackborn 1997, A = 1/lambda, B=1/lambda (still wall at x=-1, full magnitude at x=1)
% In Weldon 2008, A = 1, B=0 (uniform background flow)
%
% When B = 0, the coordinate transformation
% x -> x, y -> y - K sin(\lambda t)
% transforms to a comoving frame of the rotor. Constant
% K = A*epsilon
% can be evaluated as Kcom.
%
% In summary, the non-dimensional parameters that are
% available are:
%
% epsilon (amplitude of periodic perturbation)
% lambda  (frequency of periodic perturbation)
% c       (the x in [-1,1] coordinate of the rotor)
% tau     (period of periodic perturbation -- 2pi/lambda)
%
% Hackborn (1997) lists how these parameters correspond to physical
% parameters:
% epsilon = V * h/(2 a^2 omega)
% lambda  = h^2 alpha / a^2 omega
% c       = C / h
%
% Where the physical parameters are
% h       - half the width of the channel
% V       - amplitude of the velocity of the plate oscillation
% d = V/alpha - amplitude of the plate oscillation
% a       - rotor radius
% omega   - angular velocity of the rotor
% alpha   - angular frequency of the plate oscillation
% C       - distance of the rotor from the centerline of the channel
%
%


  properties

    epsilon % (nondimensional) amplitude of wall oscillation

    lambda  % angular frequency of wall oscillation

    c       % rotor location (between -1 and 1)

    A       % A+Bx is the profile of the oscillating velocity flow
    B       % A+Bx is the profile of the oscillating velocity flow
  end

  properties (Dependent = true, SetAccess = private)
    tau    % compute/set period of the wall oscillation
    Kcom   % comoving constant
  end


  properties (SetAccess = immutable)

    quadk   % coordinate points

    quadw   % weights

  end

  methods

    function K = get.Kcom(obj)
    % COMOVING CONSTANT

      if (obj.B ~= 0 )
        warning('Comoving frame is not valid as B is not zero.')
      end

      K = obj.A*obj.epsilon;

    end

    function out = get.tau(obj)
      out = 2*pi/obj.lambda;
    end

    function lambda = set.tau(obj, p)
      lambda = 2*pi/p;
      obj.lambda = lambda;
    end

    function obj = HackbornWeldonRotOsc( dt, flowp )
    %HACKBORNWELDONROTOSC Construct a Hackborn Rotor-Oscillator flow
    % HackbornWeldonRotOsc( dt, params )
    %
    % dt    time discretization step
    % flowp
    %     -- 1 x 3 vector of coefficients [epsilon, lambda, c, a, b]
    %                  epsilon -- strength of wall oscillation
    %                  lambda  -- angular frequency of wall oscillation
    %                  c       -- rotor location (between -1 and 1)
    %                  a,b     -- linear background cross-channel velocity
    %                             profile ( a + b x )
    %     -- 'regular'      - parameter set [0.04, 2.463, 0.54, 1/lambda, 1/lambda]
    %     -- 'structured'   - parameter set [0.02, 1.232, 0.54, 1/lambda, 1/lambda]
    %     -- 'mixing'       - parameter set [0.02, 0.406, 0.54, 1/lambda, 1/lambda]
    %     -- 'margaux'      - [0.125, 0.4*pi, 0.54, 1, 0]
    %

      if nargin < 2
        help ContinuousFlows.HackbornWeldonRotOsc.HackbornWeldonRotOsc
      end

      obj.Domain = [-2,2; -1,1];
      obj.dt = dt;

      if ischar( flowp )
        switch flowp
          case 'regular'
            flowp = [0.04, 2.463, 0.54, 1/2.463, 1/2.463];
          case 'structured'
            flowp = [0.02, 1.232, 0.54, 1/1.232, 1/1.232];
          case 'mixing'
            flowp = [0.1, 0.406, 0.54, 1/0.406, 1/0.406];
          case 'margaux'
            %            flowp = [0.125, 0.4*pi, 0.54, 1, 0];
            flowp = [0.125, 0.4*pi, 0.54, 1, 0];
          otherwise
            error('Unknown parameter set');
        end
      end

      flowp = num2cell(flowp);
      [obj.epsilon, obj.lambda, obj.c, obj.A, obj.B] = deal(flowp{:});

      %% Gauss-Legendre points and weights
      % on the k = [0,50] interval

      N = 50; % truncation order

      % quadk - column vector
      % quadw - row vector
      [K,W] = legendreweights(N, 0, 50);
      obj.quadk = K;
      obj.quadw = W;

      %% Set up integration parameters
      obj.integrator = @ode23t;
      obj.intprops = odeset;
      obj.intprops = odeset(obj.intprops, 'Vectorized', 'on');
      obj.intprops = odeset(obj.intprops, 'Jacobian', @obj.jacobian);
      obj.intprops = odeset(obj.intprops, 'MaxStep', obj.tau*5/100);
      %obj.intprops = odeset(obj.intprops, 'Stats','on' );

    end

    function [out] = Psi( obj, t, x, order )
    %% Time-varying stream function
    %  Change order to return either value, first, or second derivatives
    %
    %  second derivatives are sorted as
    %  [xx; xy; yy]

    % oscillate horizontal variable
    x(1,:,:) = x(1,:,:) - obj.epsilon*sin( obj.lambda * t );

    % swap X and Y variables to have the domain horizontal
    out = -flipud(obj.Phi(x,order) + obj.Gamma(x,order));

    end

    function [out] = Phi( obj, x, order )
    %% Closed-form (logarithmic) static term in stream function
    %  Change order to return either value, first, or second derivatives
    %
    %  second derivatives are sorted as
    %  [xx; xy; yy]

      C = obj.c;

      % swap X and Y variables to have the domain horizontal
      X = x(2,:);
      Y = x(1,:);

      CpiY = cosh(pi*Y/2);
      Cneg = cos(pi*(C-X)/2) - CpiY;
      Cpos = cos(pi*(C+X)/2) + CpiY;

      if order == 0
        out = log( -Cneg./Cpos )/2;
      elseif order == 1
        Den = Cneg.*Cpos;
        dFX =  (pi/4).*( sin(C.*pi)-2.*cos(C.*pi/2).*cosh(pi.*Y/2).*sin(pi.*X/2) )./Den;
        dFY = -(pi/2).*( cos(C.*pi/2).*cos(pi.*X/2).*sinh(pi.*Y/2) )./Den;

        sel = abs(X-C) < 1e-7 & abs(Y) < 1e-7;
        if any(sel)
          dFX(sel) = (pi/4) * tan(C*pi/2);
          dFY(sel) = 0;
        end
        out = [dFX; dFY];
        %        warning('Need to implement test for consistency of Jacobian')
        %        warning('Need to implement abstract 2d Hamiltonian')
      elseif order == 2
        Den = 8.*(Cneg.*Cpos).^2/pi^2;
        SinX2 = sin(pi.*X/2);
        CosX = cos(pi.*X);
        CoshY = cosh(pi.*Y);

        dFXX = cos(C.*pi/2).* cos(pi.*X/2) .* ...
               (cosh(pi.*Y/2).*(-3+cos(C.*pi) + CosX + CoshY) + 4.*sin(C.*pi/2).*SinX2)...
               ./ Den;
        dFYY = -dFXX;

        dFXY = ( cos(C.*pi/2).*SinX2.*sinh(pi*Y/2).*(-3+cos(C.*pi) - CosX - CoshY) + ...
                 sin(C.*pi).*sinh(pi*Y) ) ./ Den;

        out = [dFXX;
               dFXY; ...
               dFYY ];

      else
        error('Orders > 2 not implemented');
      end

    end

    function [out] = Gamma( obj, x, order )
    %% Integral static term in stream function
    %  Change order to return either value, first, or second derivatives
    %
    %  second derivatives are sorted as
    %  [xx; xy; yy]

      C = obj.c;
      Nx = size(x,2);

      persistent TanhK CothK K K2 Pp Pn

      if isempty(Pp) || isempty(Pn)

        K = obj.quadk;
        K2 = obj.quadk.^2;

        TanhK = tanh(K);
        CothK = coth(K);

        CoshCK = cosh(C*K);
        SinhCK = sinh(C*K);
        Sinh2K = sinh(2*K);

        Pp = 2*CoshCK./(Sinh2K + 2*K);
        Pn = 2*SinhCK./(Sinh2K - 2*K);

      end

      out = zeros(order+1, Nx);
      for n = 1:Nx % Nx is typically 1

        % swap X and Y variables to have the domain horizontal
        X = x(2,n);
        Y = x(1,n);

        CoshKX = cosh(K*X);
        SinhKX = sinh(K*X);

        % zeroth-order terms
        G = Pp .* (TanhK .* CoshKX - X.*SinhKX) + Pn.*(CothK.*SinhKX - X.*CoshKX);
        KY = K.*Y;
        CosKY = cos(KY);

        if order == 0
          out(n) = obj.quadw * (G.*CosKY);
          continue;
        end

        SinKY = sin(KY);

        % first-order terms
        Gx = CoshKX.*(-Pn - K.*(Pp.*X - Pn.*CothK)) + ...
             SinhKX.*(-Pp - K.*(Pn.*X - Pp.*TanhK));

        if order == 1
          % Gauss-Legendre integral as an inner product with weight
          % row-vector

          %          out(1,n) = obj.quadw * (Gx.* CosKY);
          %          out(2,n) = -obj.quadw * (G .* K.*SinKY);
          out(:,n) = obj.quadw * [(Gx.* CosKY), -(G .* K.*SinKY)];
          continue;
        end

        % order == 2
        if order == 2

          Gxx = -2*K.*(Pp.*CoshKX + Pn.*SinhKX) + ...
                K2.* (-Pn.*X.*CoshKX - Pp.*X.*SinhKX +...
                      Pn.*CothK.*SinhKX + Pp.*CoshKX.*TanhK);

          %Gamma_xx
          out(1,n) = obj.quadw * (Gxx.* CosKY);
          %Gamma_xy
          out(2,n) = -obj.quadw * (Gx .* K .*SinKY);
          %Gamma_yy
          out(3,n) = -obj.quadw * (G.* K2 .*CosKY);
          continue;
        end

        error('Higher orders not implemented');
      end

    end

  end

end

function [x,w]=legendreweights(N,a,b)
% LEGENDREWEIGHTS Compute Legendre-Gauss points and weights.
%
% [x,w]=legendreweights(N,a,b)
%
% This script is for computing definite integrals using Legendre-Gauss
% Quadrature. Computes the Legendre-Gauss nodes and weights  on an interval
% [a,b] with truncation order N
%
% Suppose you have a continuous function f(x) which is defined on [a,b]
% which you can evaluate at any x in [a,b]. Simply evaluate it at all of
% the values contained in the x vector to obtain a vector f. Then compute
% the definite integral using sum(f.*w);
%
% Written by Greg von Winckel (named: lgwt) - 02/25/2004
  N=N-1;
  N1=N+1; N2=N+2;

  xu=linspace(-1,1,N1)';

  % Initial guess
  y=cos((2*(0:N)'+1)*pi/(2*N+2))+(0.27/N1)*sin(pi*xu*N/N2);

  % Legendre-Gauss Vandermonde Matrix
  L=zeros(N1,N2);

  % Derivative of LGVM
  Lp=zeros(N1,N2);

  % Compute the zeros of the N+1 Legendre Polynomial
  % using the recursion relation and the Newton-Raphson method

  y0=2;

  % Iterate until new points are uniformly within epsilon of old points
  while max(abs(y-y0))>eps


    L(:,1)=1;
    Lp(:,1)=0;

    L(:,2)=y;
    Lp(:,2)=1;

    for k=2:N1
      L(:,k+1)=( (2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1) )/k;
    end

    Lp=(N2)*( L(:,N1)-y.*L(:,N2) )./(1-y.^2);

    y0=y;
    y=y0-L(:,N2)./Lp;

  end

  % Linear map from[-1,1] to [a,b]
  x=(a*(1-y)+b*(1+y))/2;
  x = x(:); % column vector

  % Compute the weights -- row vector
  w=(b-a)./((1-y.^2).*Lp.^2)*(N2/N1)^2;
  w=w(:).';
end
