function bin_search(num::Float64, A::Array{Float64,1})
    index = searchsortedfirst(A, num)
    if index > length(A) || A[index] != num
        index = 0
    end
    return convert(Int64,index)
end