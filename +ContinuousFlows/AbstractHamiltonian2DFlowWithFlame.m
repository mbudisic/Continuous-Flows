classdef (Abstract) AbstractHamiltonian2DFlowWithFlame < ContinuousFlows.AbstractODEFlow
    %ABSTRACTHAMILTONIAN2DFLOW Abstract class specifying interface of a "continuous-time" Hamiltonian system in terms of its stream function Psi with magnitude of front propagation
    
    properties
        v0 % magnitude of flame propagation
        delta % tolerance region in in which the flame is tapered to zero
    end
    
    methods (Abstract)
        
        [out] = Psi( obj, t, x, o )
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
        
    end
    
    methods
        
        function [ out ] = bump(obj, x,y)
        % BUMP Compute the bump function that is zero outside of the domain
        %      and 1 inside the domain, with delta-width transitionzone
        bumpx = ContinuousFlows.bump( x, obj.Domain(1,:), obj.delta );
        bumpy = ContinuousFlows.bump( y, obj.Domain(2,:), obj.delta );
        out = bumpx.*bumpy;
        end
        
        function [ f ] = vf( obj, t, x )
        % VF Compute the velocity field along a single
        % trajectory given by (t, x)
        % [ f ] = vf( obj, t, x )
        %
        % t   - row-vector of times
        % x   - trajectory
        %     - columns correspond to time steps
        %     - rows correspond to states
        %
        % Returns:
        % f   - evaluation of the velocity field
        %     - each f(:,i) is a dim x 1 velocity field evaluation
        %     - of the velocity field at [ t(i), x(i,:) ] point
        
        % system is Hamiltonian (has a stream function) but the flame front modifies
        % the velocity field
        
        
        % evaluate stream function and its derivatives
        
        dPsi = obj.Psi(t,x,1);
        ddPsi = obj.Psi(t,x,2);
        Psix = dPsi(1,:);
        Psiy = dPsi(2,:);
        Psixx = ddPsi(1,:);
        Psixy = ddPsi(2,:);
        Psiyy = ddPsi(3,:);
        
        % evaluate the front speed, tapered to zero at edges
        v = obj.bump(x(1,:),x(2,:))*obj.v0;
        
        % construct the velocity field after Mitchel, Mahoney (2012) Chaos 22, eqns. 3
        th = x(3,:);
        ux = Psiy + v.*sin(th);
        uy = -Psix - v.*cos(th);
        uth = -Psixy.*sin(2*th) - Psiyy.*sin(th).^2 - Psixx.*cos(th).^2;
        
        f = [ux;uy;uth];
        
        end
        
        function [ err ] = testPsi( obj, t, x, o, delta )
        %TESTPSI Compute difference between numeric and analytic derivatives of
        %        Psi.
        %
        % The intended use is to verify correctness of analytic expressions in the
        % derivatives of the stream function.
        %
        % err = obj.testPsi( t, x, o )
        % Compute the difference between Psi(t,x,o) and numerical
        % central-difference of Psi(t,x,o-1) at a single space-time point
        % (t,x).
        %
        % err = obj.testPsi( ..., delta )
        % Uses explicit spatial step delta (default is 1e-6)
        %
        % The difference is returned as err.
        
        if nargin < 5
            delta = 1e-6;
        end
        Nx = size(x,2);
        
        assert( o >= 1, 'Order has to be >= 1');
        assert( Nx == 1, 'Single point x has to be provided');
        assert( numel(t) == 1, 'Single point t has to be provided');
        
        %% analytic derivative
        aPsiD = obj.Psi( t, x, o );
        
        %% central difference
        stencil = eye(2)*delta;
        xi = [bsxfun(@plus, x, stencil), ...
            bsxfun(@minus, x, stencil) ];
        
        ti = repmat( t, [1, size(xi,2)] );
        nPsi = obj.Psi( ti, xi, o-1 );
        nPsiD = ( nPsi( :, 1:2) - nPsi( :, 3:4) )/(2*delta);
        
        if size(nPsiD,1) == 1
            nPsiD = nPsiD(:);
        elseif size(nPsiD,1) == 2
            % the off-diagonal elements can be different
            nPsiD1 = [nPsiD(1,1); nPsiD(2,1); nPsiD(2,2)];
            nPsiD2 = [nPsiD(1,1); nPsiD(1,2); nPsiD(2,2)];
            nPsiD = [nPsiD1,nPsiD2];
        end
        
        %% compute the error
        err = bsxfun( @minus, aPsiD, nPsiD );
        end
        
        
    end
end
