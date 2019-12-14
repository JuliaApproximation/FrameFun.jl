
using PGFPlotsX
import PGFPlotsX: Plot, PlotInc, Options, Axis
import BasisFunctions: _Plot, _PlotInc

for F in (:DictFun,)
    @eval begin
        Axis(F::$F, trailing...; opts...) =
            @pgf Axis({ymode="log"}, Plot(F, trailing...; opts...))
    end
end

for plot in (:Plot, :PlotInc)
    _plot = Meta.parse("_"*string(plot))
    for EXP in (:DictFun,)
        @eval begin
            $(plot)(F::$EXP; opts...) =
                $(plot)(Options(), F; opts...)
            $(plot)(f::Function, F::$EXP; opts...) =
                $(plot)(Options(), f, F; opts...)
            $(plot)(F::$EXP, f::Function; opts...) =
                $(plot)(Options(), f, F)
            $(plot)(F::$EXP, f::$EXP; opts...) =
                $(plot)(Options(), f, F; opts...)
        end
    end
    @eval begin

        function $(plot)(options::Options, F::DictFun; plot_extension=false, opts...)
            $(plot)(options, plot_extension ? Expansion(basis(dictionary(F)),coefficients(F)) : expansion(F); opts...)
        end

        $(plot)(options::Options, F::DictFun, f::Function; opts...) =
            $(_plot)(options, f, F; opts...)
        $(plot)(options::Options, f::Function, F::DictFun; opts...) =
            $(_plot)(options, f, F; opts...)
        $(plot)(options::Options, F::DictFun, f::DictFun; opts...) =
            $(_plot)(options, f, F; opts...)
        function $(_plot)(options::Options, f::Function, F::DictFun; plot_extension=false, opts...)
            if plot_extension
                $(_plot)(options, f, Expansion(basis(F),coefficients(F)); opts...)
            else
                $(_plot)(options, f, expansion(F); opts...)
            end
        end

        function $(_plot)(options::Options, f::DictFun, F::DictFun; plot_extension=false, opts...)
            if plot_extension
                $(_plot)(options, Expansion(basis(f),coefficients(f)), Expansion(basis(F),coefficients(F)); opts...)
            else
                $(_plot)(options, expansion(f), expansion(F); opts...)
            end
        end
    end
end
