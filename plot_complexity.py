import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

val_time=np.load('val_time.npy')
pol_time=np.load('pol_time.npy')

fig = plt.figure()
ax=fig.add_subplot(121,projection='3d')
ax2=fig.add_subplot(122,projection='3d')
y=np.arange(0,11,1)
x=np.arange(0,31,1)

X,Y = np.meshgrid(x,y)
print X.shape,Y.shape,val_time.shape

ax.plot_wireframe(X,Y,pol_time,rstride=1, cstride=1,color=u'b')

ax.set_xlabel('N')
ax.set_ylabel('M')
ax.set_zlabel('Time (sec)')

ax2.plot_wireframe(X,Y,val_time,rstride=1, cstride=1,color=u'b')
plt.show()
#plt.savefig('complexity.png')
#plt.show()

# fig=plt.figure()
# plt.plot(y,val_time[:,5],'r')
# plt.plot(y,pol_time[:,5],'b-')
# plt.plot(y,val_time[:,10],'r')
# plt.plot(y,pol_time[:,10],'b-')
# plt.show()
#print val_time,pol_time