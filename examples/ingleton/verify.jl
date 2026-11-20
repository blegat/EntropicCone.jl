module IngletonCertificate

# Each signed integer represents ±H(S), where the bits of |S| select variables.
# H(∅) = 0. Sorting gives a common representation for comparing rows.
_row(terms...) = Tuple(sort!([s for s in terms if !iszero(s)]))
_bit(i) = 1 << (i - 1)
_subsets(z) = (s for s in 0:z if s & z == s)

function shannon_rows()
    rows = Set{Tuple{Vararg{Int}}}()
    for n in (4, 10)
        full = (1 << n) - 1
        for i in 1:n
            x = _bit(i)
            # H(x | all other variables) ≥ 0.
            push!(rows, _row(full, -(full ⊻ x)))
            for j in (i + 1):n
                y = _bit(j)
                for z in _subsets(full ⊻ x ⊻ y)
                    # I(x; y | z) ≥ 0.
                    push!(rows, _row(x | z, y | z, -z, -(x | y | z)))
                end
            end
        end
    end
    return rows
end

function equality_rows()
    rows = Set{Tuple{Vararg{Int}}}()
    # Averaging the root entropy vector over a↔b and c↔d preserves the
    # almost-entropic cone, the Ingleton objective, and joint entropy.
    for (x, y) in ((_bit(1), _bit(2)), (_bit(3), _bit(4)))
        for z in _subsets(15 ⊻ x ⊻ y)
            push!(rows, _row(x | z, -(y | z)))
        end
    end
    # Chang (2026), arXiv:2610.05168, equation (20). Indices 1,…,10 denote
    # a,b,c,d,r,s,t,u,v,w. Each step gives the source tuple and conditioning
    # variables. Copies are appended in order, starting with r at index 5.
    steps = (
        ((2, 4), (1, 3)),
        ((3, 4), (1, 2, 6)),
        ((2,), (1, 3, 5, 6, 7)),
        ((2,), (3, 5, 6, 7, 9)),
    )
    n = 4
    for (source, conditioning) in steps
        old = (1 << n) - 1
        z = sum(_bit, conditioning)
        copies = (n + 1):(n + length(source))
        y = sum(_bit, copies)
        # The copies are independent of the old variables conditional on z.
        push!(rows, _row(old, y | z, -z, -(old | y)))
        for subtuple in 1:((1 << length(source)) - 1)
            chosen = [k for k in eachindex(source) if !iszero(subtuple & _bit(k))]
            src = sum(i -> _bit(source[i]), chosen)
            dst = sum(i -> _bit(copies[i]), chosen)
            for w in _subsets(z)
                push!(rows, _row(src | w, -(dst | w)))
            end
        end
        n += length(source)
    end
    return rows
end

"""
    verify()

Check the integer certificate for 124 Ing + 19 H(ABCD) ≥ 0, using the
four-copy-step construction of Chang and ordinary self-adhesivity.
Return the certified lower bound for the normalized Ingleton score.
"""
function verify()
    inequalities = shannon_rows()
    equalities = equality_rows()
    lhs = zeros(BigInt, (1 << 10) - 1)
    for line in eachline(joinpath(@__DIR__, "certificate.tsv"))
        if startswith(line, "#")
            continue
        end
        relation, weight, terms = split(line, '\t')
        coefficient = parse(BigInt, weight)
        row = _row(parse.(Int, split(terms))...)
        if relation == "ge"
            @assert row in inequalities && coefficient >= 0
        else
            @assert relation == "eq"
            @assert row in equalities || _row(map(-, row)...) in equalities
        end
        for s in row
            lhs[abs(s)] += sign(s) * coefficient
        end
    end
    # All rational dual multipliers were scaled by 372. The combination
    # must equal 372 Ing + 57 H(ABCD) in all 1,023 coordinates.
    target = zeros(Int, length(lhs))
    target[[3, 5, 9, 6, 10]] .= 372       # AB, AC, AD, BC, BD
    target[[1, 2, 12, 7, 11]] .= -372     # A, B, CD, ABC, ABD
    target[15] = 57
    @assert lhs == target
    return -57 // 372
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    println("Verified Ingleton lower bound ", IngletonCertificate.verify())
end
