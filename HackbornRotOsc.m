%%HACKBORNROTOSC
% Hackborn Rotor-Oscillator flow -- Hackborn et al., JFM, (1997)
%
% The flow evolves in the channel [-1,1] x [-inf, inf] although
% practically [-1,1] x [-4, 4] is enough for common parameters
%
% Stream function $\Psi(x,y,t)$ is given by three components:
% $$ \Psi(x,y,t) = \Phi(x,y) + \Gamma(x,y) + \epsilon \Lambda(t,x) $$
%
% $$ \Phi(x,y) = (1/2) \log \frac{\cosh(\pi y/2) - \cos[\pi(x-c)/2]}{\cosh(\pi
% y/2) + \cos[\pi(x+c)/2]} $$
%
% $$ \Gamma(x,y) = \int_0^\infty \cos(k y) G(x,k) $$
% where
% $$ G(x,k) = ... $$
%
% and
% $$\Lambda(t,x) = (x + x^2/2) \cos(\lambda t)


classdef HackbornRotOsc < ContinuousFlows.Hamiltonian2DFlow

  properties
    %% flow properties
    epsilon % strength of wall oscillation
    lambda  % frequency of wall oscillation
    c       % rotor location (between -1 and 1)

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
    % flowp is a 1 x 3 vector of flow parameters [epsilon, lambda, c]

      obj.dt = dt;

      flowp = num2cell(flowp);
      [obj.epsilon, obj.lambda, obj.c] = deal(flowp{:});

      %% Set up integration parameters
      obj.integrator = @ode23t;
      obj.intprops = odeset;
      obj.intprops = odeset(obj.intprops, 'Vectorized','on');

      %% Gauss-Legendre points and weights
      % on the k = [0,100] interval
      N = 100;
      % quadk - column vector
      % quadw - row vector
      [obj.quadk,obj.quadw] = ContinuousFlows.lgwt(N, 0, 50);

    end

    function [out] = Psi( obj, t, x, order )
    %% Time-varying stream function
    %  Change order to return either value, first, or second derivatives
    %
    %  second derivatives are sorted as
    %  [xx; xy; yy]
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
      Nx = size(x,2);

      x

      X = x(1,:);
      Y = x(2,:);

      CpiY = cosh(pi*Y/2);
      Cneg = cos(pi*(C-X)/2) - CpiY;
      Cpos = cos(pi*(C+X)/2) + CpiY;

      if order == 0
        out = log( -Cneg/Cpos )/2;
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
               / Den;
        dFYY = -dFXX;

        dFXY = ( cos(C.*pi/2).*SinX2.*sinh(pi*Y/2).*(-3+cos(C.*pi) - CosX - CoshY) + ...
               sin(C.*pi).*sinh(pi*Y) ) ./ Den;

        out = [dFXX; dFXY; dFYY];

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

      persistent CoshK SinhK TanhK CothK CoshCK SinhCK DenNeg DenPos

      if isempty(CoshK)
        CoshK = cosh(obj.quadk);
        SinhK = sinh(obj.quadk);
        TanhK = tanh(obj.quadk);
        CothK = coth(obj.quadk);
        CoshCK = cosh(C*obj.quadk);
        SinhCK = sinh(C*obj.quadk);
        DenNeg = obj.quadk - CoshK.*SinhK;
        DenPos = obj.quadk + CoshK.*SinhK;
        min(abs(DenNeg))
        min(abs(DenPos))
      end

      out = zeros(order+1, Nx);
      for n = 1:Nx

        X = x(1,n);
        Y = x(2,n);

        CoshKX = cosh(obj.quadk*X);
        SinhKX = sinh(obj.quadk*X);

        % zeroth-order terms
        g = SinhCK .* (X .* CoshKX - CothK .* SinhKX ) ./ DenNeg ...
            - CoshCK.*(X .* SinhKX - CoshKX .* TanhK) ./ DenPos;

        if order == 0
          out(n) = obj.quadw * (g.*cos(obj.quadk*Y));
          continue;
        end

        % first-order terms
        gx = SinhCK.*(CoshKX - obj.quadk.*CoshKX.*CothK)./DenNeg - ...
             CoshCK.*(SinhKX - obj.quadk.*SinhKX.*TanhK)./DenPos;

        KY = obj.quadk.*Y;
        if order == 1
          % Gauss-Legendre integral as an inner product with weight row-vector
          out(1,n) = obj.quadw * (gx.* cos(KY));
          out(2,n) = -obj.quadw * (g .* obj.quadk.*sin(KY));
          continue;
        end

        % order == 2
        if order == 2

          K2 = obj.quadk.^2;
          gxx = K2 .* ( -CothK .* SinhCK .* SinhKX ./ DenNeg + ...
                        CoshCK .* CoshKX .* TanhK ./ DenPos );

          %Gamma_xx
          out(1,n) = obj.quadw * (gxx.* cos(KY));
          %Gamma_xy
          out(2,n) = -obj.quadw * (gx .* obj.quadk.*sin(KY));
          %Gamma_yy
          out(3,n) = -obj.quadw * (g.* K2 .*cos(KY));

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

    whos

      Nx = size(x,2);
      if order == 0
        out = ( x(1,:) + x(1,:).^2 / 2 ) .* cos(obj.lambda * t );
      elseif order == 1
        out = [zeros(1,Nx); ...
               ( 1 + x(1,:) ) .* cos(obj.lambda*t ) ];
      else
        out = [zeros(2,Nx); cos(obj.lambda*t) ];
      end
    end


  end

end
