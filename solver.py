nCols = 20
nRows = 30

D = hooke(1, 2.1e9, 0.4)

nDofs = dofs.max()
nElements = nodeTopo.shape[0]

print "Number of elements = ", nElements
print "Max dofs           = ", nDofs

K = zeros( [nDofs, nDofs] )
f = zeros( [nDofs, 1] )

for i in range(24):
    
    ep = [1,0.1]
    Ke, fe = plante(ex[i,:], ey[i,:], ep, D)

    n1 = nodeTopo[i,0]
    n2 = nodeTopo[i,1]
    n3 = nodeTopo[i,2]
    edof = array( [dofs[n1,0],dofs[n1,1],dofs[n2,0],dofs[n2,1],dofs[n3,0],dofs[n3,1]], 'int32' )-1
    K[ix_((edof),(edof))] = K[ix_((edof),(edof))] + Ke

f[nDofs-1,0] = 10000

bc = ones(nDofs, 'bool')
bcDofs = arange(nDofs)
bc[0] = False
bc[1] = False
bc[int(nCols+1)*2] = False
bc[int(nCols+1)*2+1] = False
bc[2*int(nCols+1)*2] = False
bc[2*int(nCols+1)*2+1] = False
bc[3*int(nCols+1)*2] = False
bc[3*int(nCols+1)*2+1] = False

bcDofs = bcDofs[bc]
print bcDofs

fsys = f[bcDofs]

t = time()
a = linalg.solve(K[ix_((bcDofs),(bcDofs))], fsys);
print time()-t

print a

#print K

#spy(K)
#show()
#print Ke
#print fe