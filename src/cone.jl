export fulldim, EntropicCone, polymatroidcone, redundant, getinequalities, getextremerays, tight!

# Entropic Cone

abstract AbstractEntropicCone{N, T<:Real}

fulldim{N}(h::AbstractEntropicCone{N}) = N

type EntropicCone{N, T<:Real} <: AbstractEntropicCone{N, T}
  n::Int
  poly::Polyhedron{N, T}

  function EntropicCone(n::Int, p::Polyhedron{N, T})
    if ntodim(n) != N
      error("The number of variables does not match the dimension of the polyhedron")
    end
    new(n, p)
  end

  function EntropicCone(p::Polyhedron{N, T})
    new(dimton(N), p)
  end

  function EntropicCone(n::Int, A::Array{T,2}, equalities::IntSet)
    if ntodim(n) != size(A, 2)
      error("The dimension in n does not agree with the number of columns of A")
    end
    if !isempty(equalities) && last(equalities) > size(A, 1)
      error("Equalities should range from 1 to the number of rows of A")
    end
    ine = SimpleHRepresentation(-A, zeros(T, size(A, 1)), equalities)
    new(n, CDDPolyhedron{N, T}(ine))
  end

end

EntropicCone{T<:AbstractFloat}(n::Int, A::Matrix{T}) = EntropicCone{size(A, 2), Float64}(n, Matrix{Float64}(A), IntSet([]))
EntropicCone{T<:Real}(n::Int, A::Matrix{T}) = EntropicCone{size(A, 2), Rational{BigInt}}(n, Matrix{Rational{BigInt}}(A), IntSet([]))

#Base.getindex{T<:Real}(H::EntropicCone{T}, i) = DualEntropy(H.n, H.A[i,:], i in H.equalities) # FIXME

Base.copy{N, T<:Real}(h::EntropicCone{N, T}) = EntropicCone{N, T}(h.n, copy(h.poly))

function fulldim(h::AbstractEntropicCone)
  fulldim(h.poly)
end

function getinequalities(h::EntropicCone)
  removeredundantinequalities!(h.poly)
  ine = getinequalities(h.poly)
  if sum(abs(ine.b)) > 0
    error("Error: b is not zero-valued.")
  end
  [DualEntropy(ine.A[i,:], i in ine.linset) for i in 1:size(ine.A, 1)]
end

function getextremerays(h::EntropicCone)
  removeredundantgenerators!(h.poly)
  ext = SimpleVRepresentation(getgenerators(h.poly))
  if size(ext.V, 1) > 0
    error("Error: There are vertices.")
  end
  [PrimalEntropy(ext.R[i,:]) for i in 1:size(ext.R, 1)]
end

function push!{N, T, S}(H::AbstractEntropicCone{N, T}, h::AbstractDualEntropy{N, S})
  if H.n != h.n
    error("The dimension of the cone and entropy differ")
  end
  push!(H.poly, HRepresentation(h))
end
function push!{N, T, S}(H::EntropicCone{N, T}, h::Vector{DualEntropy{N, S}})
  push!(H.poly, HRepresentation(h))
end

function Base.intersect!{N}(h::AbstractEntropicCone{N}, ine::HRepresentation{N})
  if N != size(ine.A, 2)
    error("The dimension for the cone and the HRepresentation differ")
  end
  h.poly = intersect(h.poly, ine)
end
function Base.intersect!(h1::AbstractEntropicCone, h2::AbstractEntropicCone)
  if h1.n != h2.n
    error("The dimension for the cones differ")
  end
  intersect!(h1.poly, h2.poly)
end

function Base.intersect(h1::AbstractEntropicCone, h2::AbstractEntropicCone)
  if h1.n != h2.n
    error("The dimension for the cones differ")
  end
  typeof(h1)(h1.n, intersect(h1.poly, h2.poly))
end

function polymatroidcone(n::Integer, minimal = true)
  # 2^n-1           nonnegative   inequalities H(S) >= 0
  # n*2^(n-1)-n     nondecreasing inequalities H(S) >= H(T) https://oeis.org/A058877
  # n*(n+1)*2^(n-2) submodular    inequalities              https://oeis.org/A001788

  # Actually, nonnegative is not required and nondecreasing only for H([n]) and H([n] \ i)
  n_nonnegative   = minimal ? 0 : 2^n-1
  n_nondecreasing = minimal ? n : n*2^(n-1)-n
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
    if !minimal
      A[offset_nonnegative+cur_nonnegative, :] = nonnegative(n, I)
      cur_nonnegative += 1
    end
    for j = 1:n
      if !myin(j, I)
        if !minimal || card(I) == n-1
          A[offset_nondecreasing+cur_nondecreasing, :] = nondecreasing(n, I, set(j))
          cur_nondecreasing += 1
        end
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

function tight!(h::EntropicCone)
  # FIXME it doesn't work if I do not specify the type
  # The type is DualEntropy{N, Int} with N unspecified :(
  tightness = DualEntropy{Int(ntodim(h.n)), Int}[setequality(nondecreasing(h.n, setdiff(fullset(h.n), set(i)), set(i)), true) for i in 1:h.n]
  push!(h, tightness)
end
