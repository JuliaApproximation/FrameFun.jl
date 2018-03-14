# weighted_sum_frame.jl

"""
A WeightedSumFrame is the union of a finite number of copies of a single frame,
each weighted by function.
"""

const WeightedSumFrame{S,T} = MultiDict{AbstractArray{WeightedDict{S,T}},S,T}

WeightedSumFrame(weightfunctions,frame::Dictionary) = MultiDict(map(x->WeightedDict(frame,x),weightfunctions))

superdict(f::WeightedSumFrame) = superdict(element(f,1))
weightfunctions(f::WeightedSumFrame) = map(weightfunction,elements(f))
name(f::WeightedSumFrame) = "Frame consisting of " * string(length(weightfunctions(f))) * " weighted copies of " name(superdict(f))
