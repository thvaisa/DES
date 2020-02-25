import numpy as np
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import matplotlib.pyplot as plt
import copy
from scipy.spatial.transform import Rotation as R
import scipy 

print(scipy.__version__)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

rad_con = np.pi/180
Nsteps = 50
mu = 1.32712440018*10**20
M_sun = 1.989*10**30
s_time = 2456879.5
mu2 = mu/M_sun    
scale_f = 10**10
solar_f = 1.0
timeStep = 50

vec_x = np.array((1,0,0))
vec_y =  np.array((0,1,0))
vec_z =  np.array((0,0,1))


class Comet:
    a = 518060000000       #semi-major axis
    e = 0.64102       #eccentricity
    omega = 12.780*rad_con   #argument of perapsis
    ohm = 50.147*rad_con     #longiute of asceinding node
    i = 7.0405*rad_con       #inclination
    M0 = 313.71*rad_con     #mean anomaly
    t0 = 2456879.5      #at t0
    R = None
    t = None
    comet = None

#renen scharw
def position_at_t(obj, t,mu):
    #print("comet",obj.a,obj.e,obj.omega,obj.t0)
    #Step 1, solver mean anomaly
    delta_t = 86400*(t-obj.t0)      #JD change unit
    M = obj.M0 + delta_t*np.sqrt(mu/obj.a**3)
    M = np.mod(M,2*np.pi)


    #Step 2: solve kepler
    E_j = M
    for i in range(0,40):
        E_j1 = E_j-(E_j-obj.e*np.sin(E_j)-M)/(1-obj.e*np.cos(E_j))  
        E_j = E_j1

    E = E_j1
    if(obj.e>1): print("WEEEEA")
    #step 3: obtain true anomaly
    v = 2*np.arctan2(np.sqrt(1+obj.e)*np.sin(E*0.5), np.sqrt(1-obj.e)*np.cos(E*0.5))

    #step 4 get sitdance to th ecentrla body
    rc = obj.a*(1-obj.e*np.cos(E))

    #get position vector
    o = rc*np.array((np.cos(v),np.sin(v),0))

    o_v = np.sqrt(mu*obj.a)/rc*np.array((-np.sin(E),np.sqrt(1-obj.e**2)*np.cos(E),0.0))

    cw=np.cos(obj.omega)
    sw=np.sin(obj.omega)

    co=np.cos(obj.ohm)
    so=np.sin(obj.ohm)

    ci=np.cos(obj.i)
    si=np.sin(obj.i)


    #transform to cartesian
    vec_xyz = np.array(( o[0]*(cw*co-sw*ci*so)-o[1]*(sw*co+cw*ci*so),\
                       o[0]*(cw*so+sw*ci*co)+o[1]*(cw*ci*co-sw*so),\
                       o[0]*(sw*si)+o[1]*(cw*si)))

    


    vec_v = np.array(( o_v[0]*(cw*co-sw*ci*so)-o_v[1]*(sw*co+cw*ci*so),\
                       o_v[0]*(cw*so+sw*ci*co)+o_v[1]*(cw*ci*co-sw*so),\
                       o_v[0]*(sw*si)+o_v[1]*(cw*si)))


    return o,o_v,vec_xyz,vec_v,E,rc



#renen scharw
def get_periapsis_vector_xyz(obj, mu):
    M = 0
    E = 0

    #step 3: obtain true anomaly
    v = 2*np.arctan2(np.sqrt(1+obj.e)*np.sin(E*0.5), np.sqrt(1-obj.e)*np.cos(E*0.5))

    #step 4 get sitdance to th ecentrla body
    rc = obj.a*(1-obj.e*np.cos(E))

    #get position vector
    vec = rc*np.array((np.cos(v),np.sin(v),0))

    #vec_t = np.sqrt(mu*obj.a)/rc*np.array((-np.sin(E),np.sqrt(1-obj.e**2)*np.cos(E),0))

    #vec_t = np.sqrt(mu*obj.a)/rc*np.array((-np.sin(E),np.sqrt(1-obj.e**2)*np.cos(E),0))


    cw=np.cos(obj.omega)
    sw=np.sin(obj.omega)

    co=np.cos(obj.ohm)
    so=np.sin(obj.ohm)

    ci=np.cos(obj.i)
    si=np.sin(obj.i)




    #transform to cartesicometan
    vec_c = np.array(( vec[0]*(cw*co-sw*ci*so)-vec[1]*(sw*co+cw*ci*so),\
                       vec[0]*(cw*so+sw*ci*co)+vec[1]*(cw*ci*co-sw*so),\
                       vec[0]*(sw*si)+vec[1]*(cw*si)))

    return vec_c,None,rc,v


#renen scharw
def get_cartesian(obj, time, vec):

    o,o_v,vec_xyz,vec_v,E,rc = position_at_t(comet,  time, mu)
 
    #ax.scatter(vec_xyz[0]/scale_f,vec_xyz[1]/scale_f,vec_xyz[2]/scale_f,color="red")

    vec_periapsis,empty0,empty1,empty2 = get_periapsis_vector_xyz(comet, mu)

    k_x = vec_xyz/np.sqrt(np.dot(vec_xyz,vec_xyz))
        
    k_z = 0
    if(E>=np.pi):
        k_z = np.cross(k_peri,k_x)
    else:
        k_z = np.cross(k_x,k_peri)

    k_z = k_z/(np.sqrt(np.dot(k_z,k_z)))
    k_theta = np.cross(k_x,k_z)
    k_theta = k_theta/(np.sqrt(np.dot(k_theta,k_theta)))
    
    #R_transform = np.outer(vec_x,k_x)+np.outer(vec_y,k_theta)+np.outer(vec_z,k_z)
    R_transform = np.outer(k_x,vec_x)+np.outer(k_theta,vec_y)+np.outer(k_z,vec_z)

    print("r",R_transform)
    #print(R_transform)
    #cw=np.cos(obj.omega)
    #sw=np.sin(obj.omega)

    #co=np.cos(obj.ohm)
    #so=np.sin(obj.ohm)

    #ci=np.cos(obj.i)
    #si=np.sin(obj.i)

    #vec_test = vec_xyz/scale_f
    #ax.quiver(vec_test[0],vec_test[1],vec_test[2], k_r[0], k_r[1], k_r[2], length=10, normalize=True,color="red")
    
    #ax.quiver(vec_test[0],vec_test[1],vec_test[2], k_theta[0], k_theta[1], k_theta[2], length=10, normalize=True,color="magenta")
    
    #ax.quiver(vec_test[0],vec_test[1],vec_test[2], k_z[0], k_z[1], k_z[2], length=10, normalize=True,color="blue")



    #R1 = R.from_rotvec(-obj.ohm * np.array([0, 0, 1])).as_matrix()
    #R2 = R.from_rotvec(-obj.i * np.array([0, 1, 0])).as_matrix()
    #R3 = R.from_rotvec(-obj.omega * np.array([0, 0, 1])).as_matrix()

    vec = np.matmul(R_transform,vec)
    print(vec)
    #transform to cartesicometan
    #vec_xyz = np.array(( vec[0]*(cw*co-sw*ci*so)-vec[1]*(sw*co+cw*ci*so),\
    #                   vec[0]*(cw*so+sw*ci*co)+vec[1]*(cw*ci*co-sw*so),\
    #                   vec[0]*(sw*si)+vec[1]*(cw*si)))

    return vec





def get_dust_parameters(v,mu,G,M_sun,r_c):
    si = v[2]/np.sqrt(v[1]**2+v[2]**2)
    ci = v[1]/np.sqrt(v[1]**2+v[2]**2)

    a = (v[0]**2+v[1]**2+v[2]**2)/(mu*G*M_sun)-2/r_c
    b = r_c**2*(v[1]**2+v[2]**2)/(mu*G*M_sun)
    
    e_d = mu/np.abs(mu)*np.sqrt(1+a*b)

    q_d = b/(1+e_d)

    c = np.sqrt(q_d*(1+e_d)/(e_d**2*mu*G*M_sun))

    sinw_d = v[0]*c
    cosw_d = c*np.sqrt(v[1]**2+v[2]**2)-1/e_d

    return np.arctan2(si,ci),e_d,q_d,np.arctan2(sinw_d,cosw_d)
    


comet = Comet()
particles = []
comet_pos = []
for i in range(0,1):
    time = s_time+10
    o,o_v,vec_xyz,vec_v,E,rc = position_at_t(comet,  time, mu)



    k_r = -o/np.sqrt(np.dot(o,o))
    k_theta = np.cross(k_r,np.array((0,0,1)))    

 
    print(k_r)
    print(k_theta)
    print(vec_z)

    R_transform = np.outer(k_r,vec_x)+np.outer(k_theta,vec_y)+np.outer(vec_z,vec_z)





    vec_new_v = np.matmul(R_transform,o_v)

    #ax.scatter(vec_xyz[0]/scale_f,vec_xyz[1]/scale_f,vec_xyz[2]/scale_f)
    ax.text(vec_xyz[0]/scale_f, vec_xyz[1]/scale_f, vec_xyz[2]/scale_f, str(i), color='red')

    vec_periapsis,empty0,empty1,empty2 = get_periapsis_vector_xyz(comet, mu)
    k_peri = vec_periapsis/np.sqrt(np.dot(vec_periapsis,vec_periapsis))
    k_r = -vec_xyz/np.sqrt(np.dot(vec_xyz,vec_xyz))
    k_z = 0
    if(E>=np.pi):
        k_z = np.cross(k_peri,k_r)
    else:
        k_z = np.cross(k_r,k_peri)

    k_z = k_z/(np.sqrt(np.dot(k_z,k_z)))

    k_theta = np.cross(k_r,k_z)
    k_theta = k_theta/(np.sqrt(np.dot(k_theta,k_theta)))

    vec_test = vec_xyz/scale_f
    vec_polar = vec_new_v/10**4
    ax.quiver(vec_test[0],vec_test[1],vec_test[2], vec_polar[0]*k_r[0], vec_polar[0]*k_r[1],vec_polar[0]* k_r[2], length=10,color="red")

    ax.quiver(vec_test[0],vec_test[1],vec_test[2], vec_polar[1]*k_theta[0], vec_polar[1]*k_theta[1], vec_polar[1]*k_theta[2], length=10,color="green")

    ax.quiver(vec_test[0],vec_test[1],vec_test[2], vec_polar[2]*k_z[0], vec_polar[2]*k_z[1], vec_polar[2]*k_z[2], length=10,color="blue")


    print("asd",R_transform)



    vec_dust = vec_new_v+1.0*(1-0.3)**(0.333)*np.random.rand(3)*10**3

    i,e_d,q_d,omega = get_dust_parameters(vec_dust,solar_f,mu2,M_sun,rc)    



    particle = Comet()
    a = q_d/(1-e_d)
    particle.ohm = 0
    particle.a = a
    particle.e = e_d
    particle.omega = omega
    particle.i = i
    particle.comet = comet
    particle.t = time
    v = 2*np.pi-omega
    M = v-2*e_d*np.sin(v)+(3/4.0*e_d**2+1/8.0*e_d**4)*np.sin(2*v)-1/3.0*e_d**3*np.sin(3*v)+5/32.0*e_d**4*np.sin(4*v)
    particle.M0 = np.mod(M,2*np.pi) 
    particle.t0 = time
 
    particles.append(particle)



    

    '''
    if(i==0):    
        ax.scatter(vec_xyz[0]/scale_f,vec_xyz[1]/scale_f,vec_xyz[2]/scale_f,color="red")
    else:
        ax.scatter(vec_xyz[0]/scale_f,vec_xyz[1]/scale_f,vec_xyz[2]/scale_f)
    
    #vec_test = vec_xyz/10**10
    #ax.scatter(vec_test[0],vec_test[1],vec_test[2])
        
    #ax.quiver(0,0,0, k_peri[0], k_peri[1], k_peri[2], length=10, normalize=True,color="red")
    #ax.quiver(vec_test[0],vec_test[1],vec_test[2], k_r[0], k_r[1], k_r[2], length=10, normalize=True,color="red")
    
    #ax.quiver(vec_test[0],vec_test[1],vec_test[2], k_theta[0], k_theta[1], k_theta[2], length=10, normalize=True,color="magenta")
    
    #ax.quiver(vec_test[0],vec_test[1],vec_test[2], k_z[0], k_z[1], k_z[2], length=10, normalize=True,color="blue")
    

    
    #ax.quiver(vec_test[0],vec_test[1],vec_test[2], k_r[0], k_r[1], k_r[2], length=10, normalize=True,color="red")
    
    #ax.quiver(vec_test[0],vec_test[1],vec_test[2], k_theta[0], k_theta[1], k_theta[2], length=10, normalize=True,color="magenta")
    

    R_transform = np.outer(k_r,vec_x)+np.outer(k_theta,vec_y)+np.outer(k_z,vec_z)


    ax.quiver(0,0,0, k_peri[0], k_peri[1], k_peri[2], length=10, normalize=True,color="blue")
    




   
    #if(E<=np.pi):
    #    R_transform = np.outer(k_r,vec_x)+np.outer(k_theta,vec_y)+np.outer(k_z,vec_z)
    #else:
    
    

    #x:theta
    #y:r

    

    #vec_polar = np.matmul(R_transform,vec_c)
    #vec_polar = np.matmul(R_transform,vec_c)
    #print("k_z",k_z,np.dot(k_z,k_theta),np.dot(k_z,k_r))
    
    #vec_polar = np.matmul(R_transform,vec_c)/(5*10**3)
    vec_c = np.matmul(R_transform,vec_c)
    vec_polar = vec_c/10**4

    print("vec",vec_c)
    vec_test = vec_xyz/scale_f
    

    #ax.quiver(vec_test[0],vec_test[1],vec_test[2], k_r[0], k_r[1],k_r[2], length=10,color="red")
    #vec_polar = np.array((1,1,1))    

    ax.quiver(vec_test[0],vec_test[1],vec_test[2], vec_polar[0]*k_r[0], vec_polar[0]*k_r[1],vec_polar[0]* k_r[2], length=10,color="red")

    ax.quiver(vec_test[0],vec_test[1],vec_test[2], vec_polar[1]*k_theta[0], vec_polar[1]*k_theta[1], vec_polar[1]*k_theta[2], length=10,color="green")

    ax.quiver(vec_test[0],vec_test[1],vec_test[2], vec_polar[2]*k_z[0], vec_polar[2]*k_z[1], vec_polar[2]*k_z[2], length=10,color="blue")
    

    v_t = vec_polar[0]*k_r+vec_polar[1]*k_theta+vec_polar[2]*k_z

   

    #ax.quiver(vec_test[0],vec_test[1],vec_test[2], v_t[0], v_t[1],v_t[2], length=10,color="magenta")


    #ax.quiver(vec_test[0],vec_test[1],vec_test[2], vec_polar[1]*k_theta[0], vec_polar[1]*k_theta[1], vec_polar[1]*k_theta[2], length=10,color="green")
    #ax.quiver(vec_test[0],vec_test[1],vec_test[2], vec_polar[2]*k_z[0], vec_polar[2]*k_z[1], vec_polar[2]*k_z[2], length=10,color="blue")


    #print("????")
    #print(R_transform)
    #print(np.sqrt(np.dot(vec_polar,vec_polar)),np.sqrt(np.dot(vec_c,vec_c)),E)

    
    #plt.savefig("{:03d}.png".format(i))
    #plt.cla()
    #print(i)
    '''
#plt.show()



locations = []
for i in range(0,1):
    time = s_time+10 
    o,o_v,vec_xyz,vec_v,E,rc = position_at_t(particles[0],  time, mu*solar_f)  
    print("rc",rc)
    vec_xyz = get_cartesian(particles[0].comet,particles[0].t, vec_xyz)
    print(vec_xyz)
    locations.append(vec_xyz/(scale_f))




'''
i = 0
for xyz in locations:   
#    print(xyz[0],xyz[1],xyz[2])  
    
    c = np.array((1,0,0))*((i)/len(locations))
    ax.scatter(xyz[0],xyz[1],xyz[2],color=c)
    i=i+1



ax.scatter([-100,100],[-100,100],[-100,100])
plt.show()
'''






