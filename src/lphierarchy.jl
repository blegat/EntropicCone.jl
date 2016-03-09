function enforceadhesivity{N, ET<:Real}(h::EntropicCone{N, ET})
  lift = h
  last = 1
  for S in 0x1:UInt(N)
    for T in 0x1:UInt(N)
      if setdiff(S, T) == S # no intersection
        I = setdiff(setdiff(ntodim(h.n), S), T)
        cur = copy(h)
        push!(cur, submodulareq(cur.n, S, T, I))
        lift *= cur
        last += 1
        equalonsubsetsof!(G, 1, last, union(S, I))
        equalonsubsetsof!(G, 1, last, union(T, I))
      end
    end
  end
  unlift(lift)
end

export enforceadhesivity
