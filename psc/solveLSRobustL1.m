%
%  min  1/2 \|O + A x - y\|^2 + 1/2 rsL2 * \|x\|_2^2 + lambda * \|x\|_1 +
%  lambda2*|O|_21
%
% @author Jonathan Pokrass
function [x, funVal, ValueL, o] = solveLSRobustL1(A, y, lambda, lambda2, varargin)

  %calc solution size
  [m, n] = size(A);
  solm = n;
  soln = size(y, 2);

  defaultOpt = [];
  defaultOpt.rsL2 = 0;
  defaultOpt.x0 = zeros(solm, soln);
  defaultOpt.o0 = zeros(size(y));
  defaultOpt.tol = 0;
  defaultOpt.maxIter = 1e4;
  opts = parseOpt(defaultOpt, varargin{:});

  opts.x0 = [opts.x0;opts.o0];

  %% Starting point initialization
  A = [A, eye(m, m)];

  AT = A';
  ATy = A' * y;

  x = opts.x0;
  Ax = A * x;

  isSmallChange = 0; % this flag tests whether the gradient step only changes a little

  rsL2 = opts.rsL2;
  L = 1 + rsL2;
  % We assume that the maximum eigenvalue of A'A is over 1

  % assign xp with x, and Axp with Ax
  xp = x;
  Axp = Ax;
  xxp = zeros(solm + m, soln);

  % alphap and alpha are used for computing the weight in forming search point
  alphap = 0;
  alpha = 1;


 L0 = max(eigs(AT*A));
 L = L0;
  useLineSearch = true;
  Axy = Ax - y;
  funVal(1) = evalFun();

  for iterStep=2:opts.maxIter
    
    %% compute search point s based on xp and x (with beta)
    beta = (alphap - 1) / alpha;

    s = x + beta * xxp;

    %% line search for L and compute the new approximate solution x

    % compute the gradient (g) at s
    As = Ax + beta * (Ax-Axp);
    ATAs = AT * As;
    g = ATAs - ATy + rsL2 * s;

    xp = x;
    Axp = Ax;


    while (1)
      % let s walk in a step in the antigradient of s to get v
      % and then do the l1-norm regularized projection
      v = s - g / L;

      % L1-norm regularized projection
      x(1:solm,:) = sign(v(1:solm,:)) .* max(abs(v(1:solm,:)) - lambda / L, 0);

      % L21-norm proejction
      % 
      ov = v(solm+1:end, :);
      ovL2 = sqrt(sum(ov .^ 2, 2));
      %Perform projection and scale the row
      x(solm+1:end,:) = bsxfun(@times, max(1 - lambda2 ./(L*ovL2), 0), ov);
       
      Ax = A * x;

      v = s - x;  % the difference between the new approximate solution x
                  % and the search point s

      if useLineSearch
        Av = Ax - As;
        r_sum = v(:)' * v(:);
        l_sum = Av(:)' * Av(:);

        if (r_sum <= 1e-20)
          % no change after the gradient step
          isSmallChange = 1;
          break;
        end

        % the condition is ||Av||_2^2 <= (L - rsL2) * ||v||_2^2
        if(l_sum <= r_sum * (L - rsL2))
          break;
        else
          L = max(1.5 * L, l_sum / r_sum + rsL2);
          % fprintf('\n L=%5.6f',L);
        end
      else
        break;
      end
    end

    ValueL(iterStep) = L;

    
    %% update alpha and alphap, and check whether converge
    alphap = alpha;
    alpha= (1 + sqrt(4*alpha*alpha +1)) / 2;

    xxp = x - xp;
    Axy = Ax - y;
    funVal(iterStep) = evalFun();

    if (isSmallChange)
      break;
    end

    if (iterStep > 2)
      if (abs( funVal(iterStep) - funVal(iterStep-1) ) <= opts.tol)
        break;
      end
    end

    %% Adaptive restart
    if sum(diag((v)'*xxp)) > 0 %|| funVal(iterStep) > funVal(iterStep - 1)
      alphap = 0;
      alpha = 1;
    end
  end

    %ValueL
%    L
%    iterStep
%    funVal(iterStep)

    o = x(solm + 1: end,:);
    x = x(1:solm,:);

  function [value] = evalFun()
      xsol = x(1:solm,:);
      o = x(1+solm:end,:);
      value = Axy(:)'* Axy(:)/2 + rsL2/2 * xsol(:)'* xsol(:) + sum(abs(xsol(:) .* lambda(:)));
      for i=1:m
        value = value + lambda2 * norm(o(i,:), 2);
      end
  end
end
