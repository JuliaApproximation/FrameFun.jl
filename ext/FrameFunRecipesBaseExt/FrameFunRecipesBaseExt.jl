module FrameFunRecipesBaseExt

using FrameFun, RecipesBase
using FrameFun: AbstractFun

@recipe f(::Type{<:AbstractFun}, F::AbstractFun) = expansion(F)

end
