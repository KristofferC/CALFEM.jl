const LTYPES = ["-", "--", ":"]
const LCOLORS = ["k", "b", "m", "r"]
const LMARKS = ["o", "*", ""]

"""
    eldraw2(ex, ey, [plotpar=[1,1,0], elnum=zeros(0)])

Draws the 2D mesh defined by `ex`, `ey`.
"""
function eldraw2(ex::VecOrMat, ey::VecOrMat,
                plotpar = [1, 1, 0], elnum=zeros(0))
    @eval import Winston
    # TODO, Make it nice for elements with curved boundaries
    error_check_plotting(plotpar)

    ltype  = Int(plotpar[1])
    lcolor = Int(plotpar[2])
    lmark  = Int(plotpar[3])
    lmark += 1

    plot_string = LTYPES[ltype] * LCOLORS[lcolor] * LMARKS[lmark]

    nx = size(ex, 1)
    ny = size(ey, 1)

    center = [sum(ex, 1)/nx; sum(ey, 1)/ny]

    # TODO can't we make this a bit nicer /JB
    npoints = size(ex , 2)

    if npoints == 2
      xs = ex
      ys = ey
    else
      xs = [ex; ex[1, :]]
      ys = [ey; ey[1, :]]
    end
    p = Winston.plot(xs, ys, plot_string)
    for el in elnum
         Winston.text(center[1, el], center[2, el], string(el))
    end

    return p
end


"""
    eldisp2(ex, ey, ed, [plotpar=[1,1,0], sfac=1.0])

Draws the displaced 2D mesh defined by `ex`, `ey` and the displacements
given in `ed`.
"""
function eldisp2(ex::VecOrMat, ey::VecOrMat, ed::VecOrMat,
                plotpar = [1, 1, 0], sfac = 1.0)
    @eval import Winston
    # TODO, Make it nice for elements with curved boundaries
    error_check_plotting(plotpar)

    ltype =  Int(plotpar[1])
    lcolor = Int(plotpar[2])
    lmark =  Int(plotpar[3])
    lmark += 1

    plot_string = LTYPES[ltype] * LCOLORS[lcolor] * LMARKS[lmark]

    # TODO can't we make this a bit nicer /JB
    npoints = size(ex, 2)
    if npoints == 2
      xs = ex + sfac * ed[1:2:end, :]
      ys = ey + sfac * ed[2:2:end, :]
    else
      xs = [ex + sfac * ed[1:2:end, :]; ex[1, :] + sfac * ed[1, :]]
      ys = [ey + sfac * ed[2:2:end, :]; ey[1, :] + sfac * ed[2, :]]
    end
    p = Winston.plot(xs, ys, plot_string)

    return p
end

function error_check_plotting(plotpars)
    (Int(plotpars[1])) in (1,2,3)   || throw(ArgumentError("linetype must be 1,2,3"))
    (Int(plotpars[2])) in (1,2,3,4) || throw(ArgumentError("linecolor must be 1,2,3,4"))
    (Int(plotpars[3])) in (1,2,0)   || throw(ArgumentError("linemark must be 1,2,0"))
end
