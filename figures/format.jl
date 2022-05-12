using PyCall
# This is required for Python version to find tex installation
# using ```julia rcParams["text.usetex"] = true```
# /!\ This is platform specific
py"""
import os
# os.environ['PATH']='/usr/local/opt/ruby/bin:/usr/local/lib/ruby/gems/3.0.0/bin:/Users/victorboussange/.gem/ruby/2.6.0/bin:/usr/local/anaconda3/bin:/usr/local/anaconda3/condabin:/Applications/Visual Studio Code.app/Contents/Resources/app/bin:/usr/local/anaconda3/bin:/Applications/Julia-1.0.app/Contents/Resources/julia/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/Caskroom/mactex/2020.0407:/opt/X11/bin'
os.environ['PATH']='/usr/local/opt/ruby/bin:/usr/local/lib/ruby/gems/3.0.0/bin:/Users/victorboussange/.gem/ruby/2.6.0/bin:/usr/local/anaconda3/bin:/usr/local/anaconda3/condabin:/Applications/Visual Studio Code.app/Contents/Resources/app/bin:/usr/local/anaconda3/bin:/Applications/Julia-1.0.app/Contents/Resources/julia/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/Library/TeX/texbin:/opt/X11/bin'
"""
using PyPlot
mplt = pyimport("mpl_toolkits.axes_grid1.inset_locator")
ticker = pyimport("matplotlib.ticker")
FormatStrFormatter = ticker.FormatStrFormatter
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] = true
rcParams["text.latex.preamble"] = ["\\usepackage{nicefrac}",
                                    ### For now skipping the following as tex version is not working
                                    # "\\usepackage{DejaVuSans}",
                                    "\\usepackage{sansmath}",
                                    "\\sansmath",
                                    "\\usepackage{mathtools}",
                                    "\\usepackage{xcolor}",
                                    "\\newcommand\\rcirc{\\textbf{II}}",
                                    "\\newcommand\\bcirc{\\textbf{I}}"
                                    ]
rcParams["font.size"] = 8
@show rcParams["font.size"]
rcParams["axes.titlesize"] = 12
@show rcParams["axes.titlesize"]
rcParams["axes.labelsize"] = 12
@show rcParams["axes.labelsize"]
rcParams["xtick.labelsize"] = 8
rcParams["ytick.labelsize"] = 8
@show rcParams["xtick.labelsize"]
rcParams["legend.fontsize"] = 10
@show rcParams["legend.fontsize"]
rcParams["figure.titlesize"] = 12
@show rcParams["figure.titlesize"]

cm2inches(x) = x / 2.54

const FIGSIZE_S = (8,8) ./ 2.54 # cm divided by inches
const FIGSIZE_M = (11,11) ./ 2.54
const FIGSIZE_L = (18,22) ./ 2.54

#################################
### other utils for plotting ###
#################################

function _scale(x)
    (x .- mean(x)) ./ std(x)
end

function plotfit(X,Y,ax,pol)
    xfit = Float64.(X); yfit = Float64.(Y)
    # p = curve_fit((t,p) ->  p[1] * exp.(-p[2] * t), xfit, yfit, [0.,0.])
    p = Polynomials.fit(xfit,yfit,pol)
    xeval = sort!(xfit)[1:end]
    ax.plot(xeval,p.(xeval),c="tab:blue")
end
