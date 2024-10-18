using ApproxOperator

function check_rank(nodes::Vector{Node},n::Int)
    ndim = (n+1)*(n+2)/2
    k = zeros(ndim,ndim)
    for (i,node) in enumerate(nodes)
        j = 0
        x = node.x
        y = node.y
        for ii in 0:n
            j += 1
            jj = n-ii
            k[i,j] = x^ii*y^jj
        end
    end
    return rank(k)
end