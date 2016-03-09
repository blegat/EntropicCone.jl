export EntropicCone, polymatroidcone, redundant, getinequalities, getextremerays

# Entropic Cone

abstract AbstractEntropicCone{N, T<:Real}

type EntropicCone{N, T<:Real} <: AbstractEntropicCone{N, T}
  n::Int
  poly::Polyhedron{N, T}

  function EntropicCone(p::Polyhedron{N, T})
    new(fulldim(p), p)
  end

  function EntropicCone(n::Int, A::Array{T,2}, equalities::IntSet)
    if ntodim(n) != size(A, 2)
      error("The dimension in n does not agree with the number of columns of A")
    end
    if !isempty(equalities) && last(equalities) > size(A, 1)
      error("Equalities should range from 1 to the number of rows of A")
    end
    ine = InequalityDescription(A, zeros(T, size(A, 1)), equalities)
    new(n, polyhedron(ine))
  end

end

EntropicCone{T<:Real}(n::Int, A::Array{T,2}) = EntropicCone{size(A, 2), T}(n, A, IntSet([]))

#Base.getindex{T<:Real}(H::EntropicCone{T}, i) = DualEntropy(H.n, H.A[i,:], i in H.equalities) # FIXME

Base.copy{T<:Real}(h::EntropicCone{size(A, 2), T}) = EntropicCone{size(A, 2), T}(h.n, copy(h.poly))

function getpoly(h::AbstractEntropicCone)
  get(h.poly)
end

function fulldim(h::AbstractEntropicCone)
  fulldim(h.poly)
end

function getinequalities(h::EntropicCone)
  removeredundantinequalities!(getpoly(h))
  ine = getinequalities(getpoly(h))
  if sum(abs(ine.b)) > 0
    error("Error: b is not zero-valued.")
  end
  [DualEntropy(ine.A[i,:], i in ine.linset) for i in size(ine.A, 1)]
end

function getextremerays(h::EntropicCone)
  removeredundantgenerators!(getpoly(h))
  ext = getgenerators(getpoly(h))
  splitvertexrays!(ext)
  if size(ext.V, 1) > 0
    error("Error: There are vertices.")
  end
  [PrimalEntropy(ext.R[i,:]) for i in size(ext.R, 1)]
end

function push!{N, T<:Real}(H::AbstractEntropicCone{N, T}, h::AbstractDualEntropy{N, T})
  if H.n != h.n
    error("The dimension of the cone and entropy differ")
  end
  push!(H.poly, InequalityDescription{T}(h))
end

function Base.intersect!(h1::AbstractEntropicCone, h2::AbstractEntropicCone)
  if h1.n != h2.n
    error("The dimension for the cones differ")
  end
  intersect!(getpoly(h1), getpoly(h2))
end

function Base.intersect(h1::AbstractEntropicCone, h2::AbstractEntropicCone)
  if h1.n != h2.n
    error("The dimension for the cones differ")
  end
  typeof(h1)(h1.n, intersect(getpoly(h1), getpoly(h2)))
end

function polymatroidcone(n::Integer)
  # 2^n-1           nonnegative   inequalities H(S) >= 0
  # n*2^(n-1)-n     nondecreasing inequalities H(S) >= H(T) https://oeis.org/A058877
  # n*(n+1)*2^(n-2) submodular    inequalities              https://oeis.org/A001788
  n_nonnegative   = 2^n-1
  n_nondecreasing = n*2^(n-1)-n
  n_submodular    = 0
  if n >= 3
    n_submodular  = (n-1)*n*2^(n-3)
  elseif n == 2
    n_submodular  = 1
  end
  offset_nonnegative   = 0
  offset_nondecreasing = n_nonnegative
  offset_submodular    = n_nonnegative + n_nondecreasing
  cur_nonnegative   = 1
  cur_nondecreasing = 1
  cur_submodular    = 1
  A = zeros(Int, n_nonnegative + n_nondecreasing + n_submodular, ntodim(n))
  for j = 1:n
    for k = (j+1):n
      A[offset_submodular+cur_submodular, :] = submodular(n, set(j), set(k), 0x0)
      cur_submodular += 1
    end
  end
  for I = 1:ntodim(n)
    A[offset_nonnegative+cur_nonnegative, :] = nonnegative(n, I)
    cur_nonnegative += 1
    for j = 1:n
      if !myin(j, I)
        A[offset_nondecreasing+cur_nondecreasing, :] = nondecreasing(n, I, set(j))
        cur_nondecreasing += 1
        for k = (j+1):n
          if !myin(k, I)
            A[offset_submodular+cur_submodular, :] = submodular(n, set(j), set(k), I)
            cur_submodular += 1
          end
        end
      end
    end
  end
  EntropicCone(n, A)
end
