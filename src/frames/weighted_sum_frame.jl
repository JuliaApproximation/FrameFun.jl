# weighted_sum_frame.jl

"""
A WeightedSumFrame is the union of a finite number of copies of a single frame,
each weighted by function.
"""

const WeightedSumFrame{S,T} = MultiDict{AbstractVector{WeightedDict{S,T}},S,T}

WeightedSumFrame(weightfunctions,frame::Dictionary) = MultiDict(map(x->WeightedDict(frame,x),weightfunctions))

function hassuperdict(f::WeightedSumFrame)
    superdict = superdict(element(f,1))
    hassuperdict = true
    for wdict in elements(f)
        if superdict(wdict) != superdict
            hassuperdict = false
        end
    end
    hassuperdict
end

function superdict(f::WeightedSumFrame)
    if hassuperdict
        return superdict(element(f,1))
    else
        error("Error: This weighted sum frame does not have a superdict")
        return
    end
end

function name(f::WeightedSumFrame)
    if hassuperdict
        return "Frame consisting of " * string(numelements(f)) * " weighted copies of " * name(superdict(f))
    else
        return "Weight sum frame with " * string(numelements(f)) * "constituent frames"
    end
end

weightfunctions(f::WeightedSumFrame) = map(weightfunction,elements(f))
