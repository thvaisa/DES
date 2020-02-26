import numpy as np
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import matplotlib.pyplot as plt
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')


scale_v = 100000

A = np.loadtxt("./build/test.dat")
scale = 1e+10

B = A[1,:]

rv= A[0,:]
rv = rv/np.sqrt(np.dot(rv,rv))*scale_v
ax.quiver(A[1,0]-B[0],A[1,1]-B[1],A[1,2]-B[2],rv[0],rv[1],rv[2],color="red")

A = A[1:,:]
r = np.array((A[0,0],A[0,1],A[0,2]))
r = -r/np.sqrt(np.dot(r,r))


r = r*scale_v
ax.quiver(A[0,0]-B[0],A[0,1]-B[1],A[0,2]-B[2],r[0],r[1],r[2],normalize=False)


ax.scatter(A[0,0]-B[0],A[0,1]-B[1],A[0,2]-B[2],color="red")
ax.scatter((A[1:,0]-B[0])*scale,(A[1:,1]-B[1])*scale,(A[1:,2]-B[2])*scale,color="blue")



#maxV = scale*np.max(np.array((np.abs(A[:-52,0]-A[0,0]),np.abs(A[:-52,1]-A[0,1]),np.abs(A[:-52,2]-A[0,2]))))
maxV = scale*np.max(np.array((np.abs(A[:,0]-A[0,0]),np.abs(A[:,1]-A[0,1]),np.abs(A[:,2]-A[0,2]))))
#ax.scatter([-maxV+A[0,0],maxV+A[0,0]],[-maxV+A[0,1],maxV+A[0,1]],[-maxV+A[0,2],maxV+A[0,2]],color="white")
ax.scatter([-maxV,maxV],[-maxV,maxV],[-maxV,maxV],color="white")
plt.show()
