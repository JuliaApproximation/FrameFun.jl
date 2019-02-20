using BasisFunctions, Test
@testset begin  "Test Measures"
@test BasisFunctions.DiracCombMeasure(EquispacedGrid(10)) isa BasisFunctions.WeightedDiracCombMeasure
@test BasisFunctions.DiracCombMeasure(EquispacedGrid(10)) isa BasisFunctions.DiracCombMeasure
@test BasisFunctions.DiracCombMeasure(EquispacedGrid(10)) isa BasisFunctions.UniformDiracCombMeasure
@test BasisFunctions.DiracCombProbabilityMeasure(EquispacedGrid(10)) isa BasisFunctions.UniformDiracCombMeasure
@test BasisFunctions.DiracCombProbabilityMeasure(EquispacedGrid(10)) isa BasisFunctions.WeightedDiracCombMeasure
@test BasisFunctions.DiracCombProbabilityMeasure(EquispacedGrid(10)) isa BasisFunctions.DiracCombProbabilityMeasure
@test BasisFunctions.DiscreteMeasure(EquispacedGrid(10),rand(10)) isa BasisFunctions.WeightedDiracCombMeasure


@test !(BasisFunctions.DiracCombMeasure(EquispacedGrid(10)) isa BasisFunctions.DiracCombProbabilityMeasure)
@test !(BasisFunctions.DiracCombProbabilityMeasure(EquispacedGrid(10)) isa BasisFunctions.DiracCombMeasure)
@test !(BasisFunctions.DiscreteMeasure(EquispacedGrid(10),rand(10)) isa BasisFunctions.DiracCombProbabilityMeasure)
@test !(BasisFunctions.DiscreteMeasure(EquispacedGrid(10),rand(10)) isa BasisFunctions.UniformDiracCombMeasure)
end
