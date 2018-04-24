function enforceadhesivity{N, ET<:Real}(h::EntropyCone{N, ET})
    lift = h
    last = 1
    for S in indset(h)
        for T in indset(h)
            if setdiff(S, T) == S # no intersection
                I = setdiff(setdiff(ntodim(h.n), S), T)
                cur = copy(h)
                push!(cur, submodulareq(cur.n, S, T, I))
                lift *= cur
                last += 1
                equalonsubsetsof!(lift, 1, last, union(S, I))
                equalonsubsetsof!(lift, 1, last, union(T, I))
            end
        end
    end
    unlift(lift)
end

export enforceadhesivity
