import numpy as np


x0 = np.array((0.1, 0.2, 0.3));
x1 = np.array((1,2,3));
y0 = np.array((0.4 , 0.5 ,0.6));
y1 = np.array((4,5,6));
z0 = np.array((0.7 , 0.8 , 0.9));
z1 = np.array((7,8,9));


R_transform = np.outer(x0,x1)+np.outer(y0,y1)+np.outer(z0,z1)
print(R_transform) 
