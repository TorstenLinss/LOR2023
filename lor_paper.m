function lor_paper()
  test(4,14,1,"euler",7);
  test(4,14,1,"crank-nicolson",7);
  test(4,14,1,"extrapolation",7);
  test(4,14,1,"bdf2",7);
  test(4,14,1,"IIIC",7);
  test(4,14,1,"SDIRK",7);
end % lor_paper

function test(ks, ke, endT, method, ref)
  % generate a single table for "method"
  fprintf(1, 'method: %s\n', method);

  % compute reference solution by dG(2), polynomial degree p
  p= 31;
  alpha2= dg2(endT, 2^ke, p);

  fprintf(1, 'M & err & est & eff & $\\eta_{\\text{init}}$ & $\\eta_f$ & $\\eta_\\text{ell}^{M,0}$ & $\\eta_{\\Psi}$ & $\\eta_{\\delta\\psi}$ \\\\\n');
  rate.err= rate.est= rate.init= rate.f= rate.Psi= rate.delta_psi= rate.ell= 0;
  for k= ks:ke
    M= 2^k;
    N= 2^(k-1); x= [-N:N]'/N;

    [u1 est(k)]= para_fem(endT,   M, x, method, ref);

    err(k)= max(abs(refine(u1,ref) - eval_poly(alpha2, p, refine(x,ref))));

    if k > ks
      rate.err=       log2(err(k-1)/err(k));
      rate.est=       log2(est(k-1).est/est(k).est);
      rate.init=      log2(est(k-1).init/est(k).init);
      rate.f=         log2(est(k-1).f/est(k).f);
      if est(k).Psi ~= 0
        rate.Psi=       log2(est(k-1).Psi/est(k).Psi);
      end
      rate.delta_psi= log2(est(k-1).delta_psi/est(k).delta_psi);
      rate.ell=       log2(est(k-1).ell/est(k).ell);
    end

    fprintf(1, '%6d & %.3e (%.2f) & %.3e (%.2f) & 1/%d & %.3e (%.2f) & %.3e (%.2f) & %.3e (%.2f) & %.3e (%.2f) & %.3e (%.2f) \\\\\n', ...
        M, err(k), rate.err, est(k).est, rate.est, round(est(k).est/err(k)), ...
        est(k).init, rate.init, est(k).f, rate.f, est(k).ell, rate.ell, ...
        est(k).Psi, rate.Psi, est(k).delta_psi, rate.delta_psi);
  end
  fprintf(1, '\n');
end % test

function [u est]= para_fem(T, K, x, method, ref)
  % T      = final time
  % K      = number of time steps
  % x      = spatial mesh
  % ref    = level of refinement
  % method = numerical method to be used: "euler", "extrapolation", "bdf2",
  %             "crank-nicolson" "SDIRK" or "IIIC"

  % assemble matrizes and boundary conditions
  [stiff, mass]= fem(x);

  % define constants
  [gamma kappa0 kappa1 kappa1p kappa2 kappa2p]= green();

  est.f= 0;
  est.Psi= 0;
  est.delta_psi= 0;
  est.ell= 0;

  % approximate initial condition by nodal interpolant
  uo= u0(x); vo= uo; wo= uo;

  % compute psi^0
  psi_o= [0; mass(:,2:end-1) \ (stiff*uo - mass*f(0,x)); 0];

  % error estimation
  % refined spatial mesh to sample data
  x_ref= refine(x, ref);
  c_ref= c(x_ref);

  if strcmp(method, "IIIC")  % Lobatto IIIC
    % initialise Runge-Kutta-Radau
    % rp = knots
    % B  = Butcher table
    rp= [0; 1];
    B=  [1/2, -1/2; 1/2, 1/2];
  endif
  
  % (i)
  est.init= kappa0 * exp(-gamma*T) * ...
               max(abs(refine(uo, ref) - u0(x_ref)));

  % (iii)
  est.ell= kappa0 * exp(-gamma*T) * ...
               apost_ell(x, c_ref, uo, f(0, x_ref) + refine(psi_o, ref), ref);

  to= 0;
  for k= 1:K
    t= k*T/K;
    % if (2*floor(k/2) ~= k)
    %   t= t-T/(3*K);
    % end
    tau= t-to;

    if strcmp(method, "extrapolation")  % extrapolation based on Euler
      % one-step Euler
      v= [0; (mass(:,2:end-1) + tau*stiff(:,2:end-1)) \ (mass*(tau*f(t,x) + vo)); 0];

      % two-step Euler
      wh= [0; (mass(:,2:end-1) + .5*tau*stiff(:,2:end-1)) \ (mass*(.5*tau*f(t-.5*tau, x) + wo)); 0];
      w=  [0; (mass(:,2:end-1) + .5*tau*stiff(:,2:end-1)) \ (mass*(.5*tau*f(t, x) + wh)); 0];

      % extrapolated approximation
      u= 2*w - v;

      psi= (4*(wh-w)+v-vo)/tau;

      vo= v; wo= w;
    elseif strcmp(method, "euler")  % backward Euler method
      u= [0; (mass(:,2:end-1) + tau*stiff(:,2:end-1)) \ (mass*(tau*f(t,x) + uo)); 0];

      psi= -(u-uo)/tau;
    elseif strcmp(method, "bdf2")  % BDF-2 method
      if k == 1
        u= [0; (mass(:,2:end-1) + tau*stiff(:,2:end-1)) \ (mass*(tau*f(t,x) + uo)); 0];
        
        psi= -(u-uo)/tau;
        Duo= (u-uo)/tau;
      else
        alpha= (2*tau+tauo)/(tau+tauo);
        beta= - tau/(tau+tauo);
        u= [0; (alpha*mass(:,2:end-1)/tau + stiff(:,2:end-1)) ...
             \ (mass*(f(t,x) + alpha*uo/tau - beta*(uo-uoo)/tauo)); 0];
             
        Du= (u-uo)/tau;
        psi= - alpha*Du - beta*Duo;
        Duo= Du;
      endif
      uoo= uo;
      tauo= tau;
    elseif strcmp(method, "crank-nicolson")  % Crank-Nicolson time stepping
      u= [0; (mass(:,2:end-1) + 0.5*tau*stiff(:,2:end-1)) ...
            \ (mass*(tau*(f(t,x)+f(t-tau,x))/2 + uo) - 0.5*tau*stiff*uo); 0];

      psi= - psi_o - 2*(u-uo)/tau;
    elseif strcmp(method, "IIIC")  % Lobatto IIIC
      w= (kron(eye(2), mass(:,2:end-1)) + kron(B, stiff(:,2:end-1))*tau) ...
           \ (kron([1; 1], mass*uo) ...
               + tau*reshape(mass*(f(t+(rp-1)*tau,x))*B', [], 1));
      v= [0; w(1:end/2); 0];
      u= [0; w(end/2+1:end); 0];
      
      psi= (v-u)/tau;    
    elseif strcmp(method, "SDIRK")  % SDIRK
      gamma= 1-1/sqrt(2);
      k1= [0; ((mass(:,2:end-1) + gamma*tau*stiff(:,2:end-1))) \ ...
            (-stiff*uo + mass*(gamma*f(t,x)+(1-gamma)*f(to,x))); 0];
      k2= [0; ((mass(:,2:end-1) + gamma*tau*stiff(:,2:end-1))) \ ...
            (-stiff*(uo+(1-2*gamma)*tau*k1) + mass*((1-gamma)*f(t,x)+gamma*f(to,x))); 0];
      
      u= uo + tau*(k1+k2)/2;
      psi= - k1 + (k1-k2)/(2*gamma);    
    else
      error("unknown method: %s\n", method);
    endif

    % compute error estimators
    sigmaj= exp(-gamma*(T-t));

    % (ii), data oscillations
    pwlsimp= max(abs(f(t,x_ref)-2*f((t+to)/2,x_ref)+f(to,x_ref)))*tau/3;
    est.f= est.f + sigmaj*kappa0*pwlsimp;

    % psi= [0; mass(:,2:end-1) \ (stiff*u-mass*f(t,x)); 0]; % definition of psi

    % (iii), ell
    est.ell= est.ell + kappa0 * sigmaj * ...
                apost_ell(x, c_ref, u-uo, f(t, x_ref) - f(to, x_ref) ...
                         + refine(psi - psi_o, ref), ref);

    % (iv), Psi
    if !strcmp(method, "crank-nicolson")
      est.Psi= est.Psi + sigmaj*tau*max(abs((psi+psi_o)/2 + (u-uo)/tau));
    endif

    % (v)
    if k < K
      sigma_j= .25*tau*tau + .5*(T-t)*(tau + (T-to)*log((T-t)/(T-to)));
      mu= min(.25*kappa0*tau*tau, kappa1*sigma_j+kappa1p*tau*tau*tau/12);
    else
      mu= .25*kappa0*tau*tau;
    end
    est.delta_psi= est.delta_psi + sigmaj*mu*max(abs(psi-psi_o))/tau;

    % updates for next step
    psi_o= psi; uo= u; to= t;
  end % for k
  est.ell= est.ell + apost_ell(x, c_ref, u, f(t, x_ref) + refine(psi, ref), ref);

  est.est= est.init + est.f + est.ell + est.Psi + est.delta_psi;
end % para_fem

function [stiff, mass] = fem(x)
  % initialise system matrizes
  % x - spatial mesh

  N= length(x)-1;   % number of mesh intervals

  pic= c(x);  % reaction coefficient at mesh points

  % precompute mesh step sizes
  hi= x(2:end-1) - x(1:end-2);
  hip1= x(3:end) - x(2:end-1);

  % stiffness matrix
  stiff= spdiags([ - 1./hi + hi.*pic(1:end-2)/6, ...
               1./hi + 1./hip1 + (hi+hip1).*pic(2:end-1)/3, ...
               - 1./hip1 + hip1.*pic(3:end)/6], ...
             [0 1 2], N-1, N+1);

  % mass matrix
  mass= spdiags([ hi/6, (hi+hip1)/3, hip1/6 ], [0 1 2], N-1, N+1);
end % fem

function est= apost_ell(x, c_ref, u, f_ref, ref)
  q= c_ref.*refine(u, ref) - f_ref;
  Iq= q(1:ref:end);
  h= x(2:end)-x(1:end-1);

  est= max(h.^2.*max(abs([Iq(2:end), Iq(1:end-1)])')')/4 ...
        + max(abs(q-refine(Iq, ref))./c_ref);
end % apost_ell

function [z_ref]= refine(z, ref)
  z_ref= zeros(ref*(length(z)-1)+1,1);

  for j= 0:ref
    z_ref(j+1:ref:end-ref+j) = ((ref-j)*z(1:end-1) + j*z(2:end))/ref;
  end % j
  z_ref(end)= z(end);
end % refine

function alpha= dg2(T, M, p)
  % dG(2) = Runge-Kutta-Radau IIA
  % T  = final time
  % M  = number of time steps
  % p  = polynomial degree of basis functions in space

  % initialise Runge-Kutta-Radau
  % rp = Radau points
  % B  = Butcher table
  rp= [2/5-sqrt(6)/10; 2/5+sqrt(6)/10; 1];
  B=  [11/45-7*sqrt(6)/360, 37/225-169*sqrt(6)/1800, -2/225+sqrt(6)/75; ...
       37/225+169*sqrt(6)/1800, 11/45+7*sqrt(6)/360, -2/225-sqrt(6)/75; ...
       4/9-sqrt(6)/36, 4/9+sqrt(6)/36, 1/9];

  % collocation at Chebychev points
  % p+1 basis functions ==> p-1 collocation points
  cp= cos((2*[1:p-1]'-1)*pi/(2*(p-1)));

  % basis functions at collocation points
  [u up upp]= basis(cp, p);

  % compute representation of initial condition
  alpha= u \ u0(cp);

  % precompute differential operator at collocation points
  L = -diag(a(cp))*upp + diag(b(cp))*up + diag(c(cp))*u;

  to= 0;
  for k= 1:M
    t= k*T/M;
    tau= t-to;

    alpha= [(kron(eye(3), u)/tau + kron(B, L)) \ ...
            (kron([1; 1; 1], u*alpha)/tau ...
               + reshape(f(t+(rp-1)*tau,cp)*B', [], 1)) ...
           ](2*end/3+1:end);
    to= t;
  end
end % dG2

function u= eval_poly(alpha, p, x)
  % basis functions at plot points
  % alpha - coefficients in basis coordinates
  % p     - polynomial degree
  % x     - points for evaluating polynomial
  [phi phip phipp]= basis(x, p);

  u= phi*alpha;
end % eval_poly

function [u up upp]= basis(x, p)
  % - basis of integrated Legendre polynomials of degree 2...p
  % - this way Dirichlet bcs are built in
  % - compute basis functions and their 1st and 2nd order
  % derivatives on the mesh x

  % Legendre polynomials and their derivatives at x
  P= zeros(length(x), p+1); Pp= P;
  P(:,1)= 1; P(:,2)= x; Pp(:,2)= 1;
  for i = 1:p-1
    P(:,i+2)= ((2*i+1)*x.*P(:,i+1) - i*P(:,i))/(i+1);
    Pp(:,i+2)= ((2*i+1)*(x.*Pp(:,i+1)+P(:,i+1)) - i*Pp(:,i))/(i+1);
  end

  % integrated Legendre polynomials and their derivatives at x
  u= zeros(length(x), p-1); up= u; upp= u;
  for i = 2:p
    u(:,i-1)= (P(:,i+1)-P(:,i-1))/(2*i-1);
    up(:,i-1)= P(:,i);
    upp(:,i-1)= Pp(:,i);
  end
end % basis

function u0= u0(x)
  % initial condition
  % x - mesh

  u0= sin(pi*((x+1)/2));
end % u0

function a= a(x)
  % diffusion coefficient
  % x mesh

  a= ones(size(x));
end % a

function b= b(x)
  % convection coefficient
  % x mesh

  b= zeros(size(x));
end % b

function c= c(x)
  % reaction coefficient
  % x mesh

  c= 5*x+6;
end % c

function f= f(t,x)
  % source term
  % t - time
  % x - mesh

  f= zeros(length(x), length(t));
  for i= 1:length(t)
    f(:,i)= exp(-4*t(i)) + cos(pi*(x+t(i)).^2);
  end
end % f

function [gamma kappa0 kappa1 kappa1p kappa2 kappa2p]= green()
  % constants in Green-function bounds

  gamma= 1/2;
  kappa0= 1;
  kappa1= 3*2^(-3/2);
  kappa1p= 0;
  kappa2= 36*kappa1;
  kappa2p= 0;
end % green
