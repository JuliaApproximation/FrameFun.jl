#######################
# Parameter sequences
#######################

"""
A `DimensionSequence` is a sequence of dimensions.
It can be indexed with an integer.
"""
abstract type DimensionSequence
end

"A doubling sequence with a given initial value."
struct DoublingSequence <: DimensionSequence
    initial ::  Int
end

# We arbitrarily choose a default initial value of 2
DoublingSequence() = DoublingSequence(2)

initial(s::DoublingSequence) = s.initial

getindex(s::DoublingSequence, idx::Int) = initial(s) * (2<<(idx-2))


"A multiply sequence with a given initial value."
struct MultiplySequence <: DimensionSequence
    initial ::  Int
    t       ::  Real
end

# We arbitrarily choose a default initial value of 2
MultiplySequence() = DoublingSequence()

initial(s::MultiplySequence) = s.initial

getindex(s::MultiplySequence, idx::Int) =  round(Int, initial(s) *s.t^(idx-1))

"A stepping with a given initial value."
struct SteppingSequence <: DimensionSequence
    initial ::  Int
    step    ::  Int
end

# We arbitrarily choose a default initial value of 2
SteppingSequence() = SteppingSequence(1,1)
SteppingSequence(init::Int) = SteppingSequence(init,1)
initial(s::SteppingSequence) = s.initial

getindex(s::SteppingSequence, idx::Int) = initial(s) + s.step*(idx-1)


"A tensor product sequences with given initial values."
struct TensorSequence
    sequences
end

getindex(s::TensorSequence, idx::Int) = [si[idx] for si in s.sequences]
