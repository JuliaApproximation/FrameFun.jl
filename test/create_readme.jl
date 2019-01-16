
println("Create README.md")

import FrameFun
jupyter = chomp(read(pipeline(`find $(homedir())/.julia/packages/Conda/  -name "jupyter" -type f`,`head -n 1`),String))
FRAMEFUNSRC = pathof(FrameFun)
FRAMEFUNPATH = splitdir(splitdir(FRAMEFUNSRC)[1])[1]
run(`$jupyter nbconvert --execute --to markdown --output $FRAMEFUNPATH/README.md $FRAMEFUNPATH/README.ipynb --ExecutePreprocessor.timeout=180`)
