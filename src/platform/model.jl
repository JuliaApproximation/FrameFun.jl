"""
A `ModelPlatform` is a platform based on a model dictionary. The platform is
defined by resizing the dictionary, using its own implementation of `resize`.
All other operations are the defaults for the model dictionary.

This platform is convenient to compute adaptive approximations based on an
example of a dictionary from the desired family.
"""
struct ModelPlatform <: Platform
    model   ::  Dictionary
end

model(p::ModelPlatform) = p.model

dictionary(p::ModelPlatform, n) = resize(model(p), n)

param_first(p::ModelPlatform) = dimensions(model(p))

SamplingStyle(p::ModelPlatform) = SamplingStyle(model(p))
SolverStyle(p::ModelPlatform, samplingstyle::SamplingStyle) = SolverStyle(model(p), samplingstyle)

measure(platform::ModelPlatform) = measure(model(platform))

elements(platform::ModelPlatform) = map(ModelPlatform, elements(model(platform)))

element(platform::ModelPlatform, i) = ModelPlatform(element(model(platform), i))
