"Return a list of all indices in an n-dimensional hyperbolic cross."
function index_set_hyperbolic_cross(s, n, α = 1)
    I = _index_set_hyperbolic_cross(s, n, α)
    [tuple((1+i)...) for i in I]
end

function _index_set_hyperbolic_cross(s, n, α = 1)
    if n == 1
        smax = floor(Int, s^(1/α))-1
        I = [[i] for i in 0:smax]
    else
        I = Array{Array{Int,1}}(0)
        I_rec = _index_set_total_degree(s, n-1)
        for idx in I_rec
            for m in 0:floor(Int,s^(1/α)/prod(1+abs.(idx)))-1
                push!(I, [idx...; m])
            end
        end
        I
    end
end
