classdef HackbornRotOsc < ContinuousFlows.AbstractHamiltonian2DFlow
%HACKBORNROTOSC Hackborn Rotor-Oscillator flow -- Hackborn et al., JFM, (1997)
%
% The flow evolves in the channel [-1,1] x [-inf, inf] although
% practically [-1,1] x [-4, 4] is enough for common parameters
%
% Stream function $\Psi(x,y,t)$ is given by three components:
% $$ \Psi(x,y,t) = \Phi(x,y) + \Gamma(x,y) + \epsilon \Lambda(t,x) $$
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
% Wall-induced periodic shear
% $$\Lambda(t,x) = (x + x^2/2) \cos(\lambda t)

  properties
    %% flow properties
    epsilon % strength of wall oscillation
    lambda  % angular frequency of wall oscillation
    c       % rotor location (between -1 and 1)
  end

  properties (SetAccess = immutable)

    %% quadrature
    quadk   % coordinate points
    quadw   % weights

  end

  methods

    function obj = HackbornRotOsc( dt, flowp )
    %HACKBORNROTOSC Construct a Hackborn Rotor-Oscillator flow
    % HackbornRotOsc( dt, params )
    %
    % dt    time discretization step
    % flowp
    %     -- 1 x 3 vector of coefficients [epsilon, lambda, c]
    %                  epsilon -- strength of wall oscillation
    %                  lambda  -- angular frequency of wall oscillation
    %                  c       -- rotor location (between -1 and 1)
    %     -- 'regular'      - parameter set [0.04, 2.463, 0.54]
    %     -- 'structured'   - parameter set [0.02, 1.232, 0.54]
    %     -- 'mixing'       - parameter set [0.02, 0.406, 0.54]
    %
    % Alternatively, set flowp to string 'Hackborn' to get a seto
    % of parameters from Hackborn 1997 paper.

      if nargin < 2
        help ContinuousFlows.HackbornRotOsc.HackbornRotOsc
      end

      obj.Domain = [-1,1; -2,2];
      obj.dt = dt;

      if ischar( flowp )
        switch flowp
          case 'regular'
            flowp = [0.04, 2.463, 0.54];
          case 'structured'
            flowp = [0.02, 1.232, 0.54];
          case 'mixing'
            flowp = [0.1, 0.406, 0.54];
          otherwise
            error('Unknown parameter set');
        end
      end

      flowp = num2cell(flowp);
      [obj.epsilon, obj.lambda, obj.c] = deal(flowp{:});

      %% Gauss-Legendre points and weights
      % on the k = [0,100] interval
      N = 150;
      % quadk - column vector
      % quadw - row vector
      [K,W] = ContinuousFlows.lgwt(N, 0, 50);
      obj.quadk = K;
      obj.quadw = W;


      %% Set up integration parameters
      obj.integrator = @ode23t;
      obj.intprops = odeset;
      obj.intprops = odeset(obj.intprops, 'Vectorized', 'on');
      obj.intprops = odeset(obj.intprops, 'Jacobian', @obj.jacobian);
      obj.intprops = odeset(obj.intprops, 'MaxStep', 1e-1);
      %obj.intprops = odeset(obj.intprops, 'Stats','on' );

      warning('Function Gamma should really be vectorized')
    end

    function [out] = Psi( obj, t, x, order )
    %% Time-varying stream function
    %  Change order to return either value, first, or second derivatives
    %
    %  second derivatives are sorted as
    %  [xx; xy; yy]


    %out = obj.Phi(x,order);
    %out = obj.Gamma(x,order);
    %out = obj.Lambda(t,x,order);

    out = obj.Phi(x,order) + ...
          obj.Gamma(x,order) + ...
          obj.epsilon * obj.Lambda(t,x,order);

    end

    function [out] = Phi( obj, x, order )
    %% Closed-form (logarithmic) static term in stream function
    %  Change order to return either value, first, or second derivatives
    %
    %  second derivatives are sorted as
    %  [xx; xy; yy]

      C = obj.c;

      X = x(1,:);
      Y = x(2,:);

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

        out = [dFXX; ...
               dFXY; ...
               dFYY];

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
      for n = 1:Nx

        X = x(1,n);
        Y = x(2,n);

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
          % Gauss-Legendre integral as an inner product with weight row-vector
          out(1,n) = obj.quadw * (Gx.* CosKY);
          out(2,n) = -obj.quadw * (G .* K.*SinKY);
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

    function [out] = Lambda( obj, t, x, order )
    %% Time-varying term in stream function
    %  Change order to return either value, first, or second derivatives
    %
    %  second derivatives are sorted as
    %  [xx; xy; yy]

      Nx = size(x,2);
      COS = cos(obj.lambda*t);
      if order == 0
        out = ( x(1,:) + x(1,:).^2 / 2 ) .* COS;
      elseif order == 1
        out = [( 1 + x(1,:) ) .* COS;
               zeros(1,Nx) ];
      else
        if numel(COS) == 1
          COS = repmat(COS, [1,Nx]);
        end
        out = [COS; zeros(2,Nx)];
      end
    end
  end
end
