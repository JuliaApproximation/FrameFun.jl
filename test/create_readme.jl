
println("Create README.md")

import FrameFun
using Conda: PYTHONDIR
jupyter = Conda.PYTHONDIR * "jupyter"
FRAMEFUNSRC = pathof(FrameFun)
FRAMEFUNPATH = splitdir(splitdir(FRAMEFUNSRC)[1])[1]
run(`$jupyter nbconvert --execute --to markdown --output $FRAMEFUNPATH/README.md $FRAMEFUNPATH/README.ipynb --ExecutePreprocessor.timeout=180`)
