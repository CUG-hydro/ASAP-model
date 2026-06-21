export find_jwt, flowdir, zbrent


"""
    zbrent(func, x1, x2; tol=1e-6, maxiter=100)

Brent's method for root finding.

Finds a root of `func(x) = 0` in the bracketing interval `[x1, x2]`, i.e. an `x`
where `sign(func(x1)) != sign(func(x2))`. The method combines inverse quadratic
interpolation, secant steps and bisection, guaranteeing convergence as long as
the root is bracketed (Numerical Recipes §9.3).

Translation of `fortran/module_initial.f90:638 zbrent`.

# Arguments
- `func`: callable `f(x)::Real` (must be cheap to evaluate; called many times)
- `x1`, `x2`: bracketing interval (signs at the endpoints must differ)
- `tol`: absolute convergence tolerance on `x` (default `1e-6`)
- `maxiter`: maximum Brent iterations (default `100`)

# Returns
- `x::Float64`: refined root location

# Throws
- `ArgumentError` if `func(x1)` and `func(x2)` have the same sign (root not bracketed)
- `ErrorException` if convergence not reached after `maxiter` iterations

# Reference
- Press, W. H. et al. (2007) Numerical Recipes: The Art of Scientific Computing,
  §9.3 Brent's Method.
"""
function zbrent(func, x1, x2; tol::Real=1e-6, maxiter::Int=100)
  a = Float64(x1)
  b = Float64(x2)
  fa = Float64(func(a))
  fb = Float64(func(b))
  EPS = eps(Float64(a))

  if (fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)
    throw(ArgumentError("zbrent: root must be bracketed (sign(func(x1)) == sign(func(x2)))"))
  end

  c = b
  fc = fb
  d = b - a
  e = d

  for _ in 1:maxiter
    if (fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)
      c = a
      fc = fa
      d = b - a
      e = d
    end
    if abs(fc) < abs(fb)
      # Sequential swap (NOT Julia's parallel tuple assignment), to match
      # the Fortran `a = b; b = c; c = a` ordering exactly:
      #   a_new = b_old, b_new = c_old, c_new = b_old (= a_new)
      tmp_a = a
      a = b
      b = c
      c = a  # = b_old
      tmp_fa = fa
      fa = fb
      fb = fc
      fc = fa  # = fb_old
      _ = (tmp_a, tmp_fa)  # silence unused warnings
    end
    tol1 = 2.0 * EPS * abs(b) + 0.5 * tol
    xm = 0.5 * (c - b)
    if abs(xm) <= tol1 || fb == 0.0
      return b
    end
    if abs(e) >= tol1 && abs(fa) > abs(fb)
      s = fb / fa
      if a == c
        p = 2.0 * xm * s
        q = 1.0 - s
      else
        q = fa / fc
        r = fb / fc
        p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0))
        q = (q - 1.0) * (r - 1.0) * (s - 1.0)
      end
      if p > 0.0
        q = -q
      end
      p = abs(p)
      if 2.0 * p < min(3.0 * xm * q - abs(tol1 * q), abs(e * q))
        e = d
        d = p / q
      else
        d = xm
        e = d
      end
    else
      d = xm
      e = d
    end
    a = b
    fa = fb
    b = b + (abs(d) > tol1 ? d : copysign(tol1, xm))
    fb = Float64(func(b))
  end

  throw(ErrorException("zbrent: exceeded maximum iterations ($maxiter)"))
end


function find_jwt(wtd::FT, z₋ₕ::Vector{FT}) where {FT<:Real}
  nzg = length(z₋ₕ) - 1
  jwt = 1
  for k in 1:nzg
    if wtd < z₋ₕ[k]
      jwt = k # 地下水水位所在的上一层
      break
    end
  end
  return jwt
end


"""
    flowdir(流向参数...)

根据流向编码确定下游网格位置

流向编码：
- 1: 东 (→)
- 2: 东南 (↘)
- 4: 南 (↓)
- 8: 西南 (↙)
- 16: 西 (←)
- 32: 西北 (↖)
- 64: 北 (↑)
- 128: 东北 (↗)
"""
function flowdir(fd::Matrix{Int}, ii::Int, jj::Int)
  # 确定j方向
  j = if fd[ii, jj] in [2, 4, 8]
    jj - 1
  elseif fd[ii, jj] in [1, 16]
    jj
  elseif fd[ii, jj] in [32, 64, 128]
    jj + 1
  else
    0
  end

  # 确定i方向
  i = if fd[ii, jj] in [128, 1, 2]
    ii + 1
  elseif fd[ii, jj] in [4, 64]
    ii
  elseif fd[ii, jj] in [8, 16, 32]
    ii - 1
  else
    0
  end
  return i, j
end
