export singleton, set, card, fullset

const EntropyIndex = UInt

# Set Manipulation
# function subsets(S)
#   if S == 0
#     return []
#   else
#     cum = [S]
#     for i = 0:(n-1)
#       if (S & (1 << i)) != 0
#         append!(cum, subsets(S $ (1 << i)))
#       end
#     end
#     return cum
#   end
# end
#
# function subsetsones(S)
#   ret = zeroentropy(N)
#   ret[subsets(S)] = 1
#   return ret
# end

function Base.setdiff(S::EntropyIndex, I::EntropyIndex)
    S & (~I)
end

# S ∪ T
function Base.union(S::EntropyIndex, T::EntropyIndex)
    S | T
end

# S ∩ T
function Base.intersect(S::EntropyIndex, T::EntropyIndex)
    S & T
end

# S ⊆ T
function issubset(S::EntropyIndex, T::EntropyIndex)
    Base.setdiff(S, T) == zero(S)
end

function singleton(i::Integer)
    EntropyIndex(1) << (i-1)
end

function fullset(n::Integer)
    (EntropyIndex(1) << n) - EntropyIndex(1)
end

function set(i::Integer)
    ret = emptyset()
    while i > 0
        ret = union(ret, singleton(i % 10))
        i = div(i, 10)
    end
    ret
end
function set(I::AbstractArray{S}) where S<:Integer
    ret = emptyset()
    for i in I
        ret = union(ret, singleton(i))
    end
    ret
end


function myin(i::Signed, I::EntropyIndex)
    (singleton(i) & I) != 0
end

function card(S::EntropyIndex)
    sum = 0
    while S > 0
        sum += S & 1
        S >>= 1
    end
    Signed(sum)
end

function emptyset()
    EntropyIndex(0)
end

function setsto(J::EntropyIndex)
    EntropyIndex(1):J
end

function mymap(map::Vector, S::EntropyIndex, n)
    T = emptyset()
    for i in 1:n
        if myin(i, S)
            T = T ∪ singleton(map[i])
        end
    end
    T
end
