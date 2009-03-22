def mlscalar2d(coords, edof, a):
    if not haveMlab:
        return

    x = reshape(coords[:,0],[size(coords[:,0])])
    y = reshape(coords[:,1],[size(coords[:,1])])
    z = zeros([size(coords[:,0])])
    ascalar = reshape(asarray(a),[size(a)])

    mlab.triangular_mesh(x, y, z, edof-1, scalars=ascalar, representation="surface")

def mlflux2d(coords, vf, scalefactor=None, displaymode="2darrow"):
    if not haveMlab:
        return

    x = reshape(coords[:,0],[size(coords[:,0])])
    y = reshape(coords[:,1],[size(coords[:,1])])
    z = zeros([size(coords[:,0])])
    u = reshape(vf[:,0],[size(vf[:,0])])
    v = reshape(vf[:,1],[size(vf[:,1])])
    w = zeros([size(vf[:,0])])

    if scalefactor == None:
        mlab.quiver3d(x, y, z, u, v, w, mode=displaymode)
    else:
        mlab.quiver3d(x, y, z, u, v, w, mode=displaymode, scale_factor=scalefactor)


def mlwireframe2d(coords, edof):
    if not haveMlab:
        return

    x = reshape(coords[:,0],[size(coords[:,0])])
    y = reshape(coords[:,1],[size(coords[:,1])])
    z = zeros([size(coords[:,0])])+1
    scalars = ones([size(coords[:,0])])

    mlab.triangular_mesh(x, y, z, edof-1, scalars=scalars, representation="mesh", colormap="bone", line_width=20.0)
    

def eldisp2(ex,ey,ed,rat=0.2):
    nen = -1
    if ex.shape != ey.shape:
        print "ex and ey shapes do not match."
        return 1.0
    
    dlmax = 0.
    edmax = 1.
    
    print rank(ex)

    if rank(ex)==1:
        nen = ex.shape[0]
        nel = 1
        dxmax=ex.T.max()-ex.T.min()
        dymax=ey.T.max()-ey.T.min()
        dlmax=max(dxmax,dymax);
        edmax=abs(ed).max();
    else:
        nen = ex.shape[1]
        nel = ex.shape[0]
        dxmax=ex.T.max()-ex.T.min()
        dymax=ey.T.max()-ey.T.min()
        dlmax=max(dxmax,dymax);
        edmax=abs(ed).max();
        
    k = rat
    return k*dlmax/edmax
