import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


from math import sin, cos
from numpy import sqrt

from matplotlib import transforms
from matplotlib.patches import Ellipse


def theta2R(theta):
    return np.array([[cos(theta), -sin(theta)],
                     [sin(theta), cos(theta)]]) 

def cov_ellipse_from_meancov(mean, cov, ax, n_std=3.0, facecolor='none', **kwargs):

    pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
    # Using a special case to obtain the eigenvalues of this
    # two-dimensionl dataset.
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2,
                      facecolor=facecolor, **kwargs)

    # Calculating the stdandard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    scale_x = np.sqrt(cov[0, 0]) * n_std
    mean_x = mean[0]

    # calculating the stdandard deviation of y ...
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_y = mean[1]

    transf = transforms.Affine2D() \
        .rotate_deg(45) \
        .scale(scale_x, scale_y) \
        .translate(mean_x, mean_y)

    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)

def cov_ellipse_from_samples(x, y, ax, n_std=3.0, facecolor='none', **kwargs):
    """
    Create a plot of the covariance confidence ellipse of *x* and *y*.

    Parameters
    ----------
    x, y : array-like, shape (n, )
        Input data.

    ax : matplotlib.axes.Axes
        The axes object to draw the ellipse into.

    n_std : float
        The number of standard deviations to determine the ellipse's radiuses.

    **kwargs
        Forwarded to `~matplotlib.patches.Ellipse`

    Returns
    -------
    matplotlib.patches.Ellipse
    """
    if x.size != y.size:
        raise ValueError("x and y must be the same size")

    cov = np.cov(x, y)
    pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
    # Using a special case to obtain the eigenvalues of this
    # two-dimensionl dataset.
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2,
                      facecolor=facecolor, **kwargs)

    # Calculating the stdandard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    scale_x = np.sqrt(cov[0, 0]) * n_std
    mean_x = np.mean(x)

    # calculating the stdandard deviation of y ...
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_y = np.mean(y)

    transf = transforms.Affine2D() \
        .rotate_deg(45) \
        .scale(scale_x, scale_y) \
        .translate(mean_x, mean_y)

    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)



'''
configuration
'''
import_motion = True


'''
constant def
'''
I2 = np.eye(2)
cross2 = np.array([[0 , -1],
                   [1 , 0 ]])


if(import_motion):
    '''
    motion import
    '''
    data = np.loadtxt("m10s.txt")

    T = np.average(data[1:]-data[:-1])
    t = data[:,0]
    N = t.size

    pos = np.transpose(data[:,1:3])
    yaw = np.transpose(data[:,3])
    vb = np.transpose(data[:,4:6])
    w = np.transpose(data[:,6])

else:
    '''
    motion generation
    '''
    T = 0.01
    t = np.arange(0,10,T,dtype=float)
    N = t.size

    vb = np.zeros((2,N))
    w = np.zeros(N)
    pos = np.zeros((2,N))
    yaw = np.zeros(N)



'''
landmark generation
'''

L = 1
lpw = np.zeros((2,L))


'''
measurement generation
'''
vel_var = 0.1
w_var = 0.1
lidar_var = 0.1

from numpy.random import default_rng
rng = default_rng(1)
vx_noise = vel_var * rng.standard_normal(N)
vy_noise = vel_var * rng.standard_normal(N)
w_noise = w_var * rng.standard_normal(N) 

vb_meas = vb + np.vstack([vx_noise, vy_noise])
w_meas = w + w_noise

lx_noise = lidar_var * rng.standard_normal((L,N))
ly_noise = lidar_var * rng.standard_normal((L,N))

lpm = np.zeros((2,L,N))

for i in range(N):
    Rwb = theta2R(yaw[i])
    twb = pos[:,i]
    for j in range(L):
        lpm[:,j,i] = np.transpose(Rwb) @ (lpw[:,j] - twb)

for j in range(L):
    lpm[0,j,:] = lpm[0,j,:] + lx_noise[j,:]
    lpm[1,j,:] = lpm[1,j,:] + ly_noise[j,:]

'''
fig1, (ax1, ax2) = plt.subplots(2,1)
ax1.plot(t, vb_meas[0,:])
ax1.set_title("vbx measurement")
ax2.plot(t, vb_meas[1,:])
ax2.set_title("vby measurement")
'''

 
'''
uncertainty propagation without fix (linear Gaussain)
'''
Q = np.diag([(T*vel_var)**2, (T*vel_var)**2, (T*w_var)**2])
R = np.diag([lidar_var**2, lidar_var**2])

# set initial condition  
xr = np.zeros((3,N))
covr = np.zeros((3,3,N))
covr[:,:,0] = 1e-10 * np.eye(3)

for k in range(1,N):
    pos0 = xr[0:2,k-1]
    yaw0 = xr[2,k-1]
    vbm0 = vb_meas[:,k-1]
    Rwb0 = theta2R(yaw0)
    w0 = w_meas[k-1]

    dpos = vbm0 * T
    pos_pred = pos0 + Rwb0 @ dpos
    yaw_pred = yaw0 + w0 * T

    PHI_pos2pos = I2
    PHI_theta2pos = cross2 @ Rwb0 @ dpos
    PHI_2theta = np.array([0, 0, 1])

    PHI_2pos = np.hstack([PHI_pos2pos, PHI_theta2pos.reshape(2,1)])

    PHI = np.vstack([PHI_2pos, PHI_2theta])


    covr_pred = PHI @ covr[:,:,k-1] @ np.transpose(PHI) + Q

    # no fix in this loop , so directly use prediction value

    xr[0:2,k] = pos_pred
    xr[2,k] = yaw_pred
    covr[:,:,k] = covr_pred

fig_pose_estimate, (ax_pos_x, ax_pos_y, ax_yaw) = plt.subplots(3,1)
ax_pos_x.plot(t,pos[0,:], 'r--', label= 'x pos gt', linewidth = 2.0)
ax_pos_x.plot(t,xr[0,:], 'b', label= 'x pos estimation', linewidth = 1.0)
ax_pos_y.plot(t,pos[1,:], 'r--', label= 'y pos gt', linewidth = 2.0)
ax_pos_y.plot(t,xr[1,:], 'b', label= 'y pos estimation', linewidth = 1.0)
ax_yaw.plot(t,yaw, 'r--', label= 'yaw gt', linewidth = 2.0)
ax_yaw.plot(t,xr[2,:], 'b', label= 'yaw estimation', linewidth = 1.0)
ax_pos_x.legend()
ax_pos_y.legend()
ax_yaw.legend()

fig_pose_err, (ax_pos_err_x, ax_pos_err_y, ax_yaw_err) = plt.subplots(3,1)
ax_pos_err_x.plot(t,pos[0,:]-xr[0,:], 'r', label='x pos error', linewidth = 1.0)
ax_pos_err_x.plot(t,3*sqrt(covr[0,0,:]), 'g--', label='+3 sigma', linewidth = 0.5)
ax_pos_err_x.plot(t,-3*sqrt(covr[0,0,:]), 'g--', label='-3 sigma', linewidth = 0.5)
ax_pos_err_y.plot(t,pos[1,:]-xr[1,:], 'r', label='y pos error', linewidth = 1.0)
ax_pos_err_y.plot(t,3*sqrt(covr[1,1,:]), 'g--', label='+3 sigma', linewidth = 0.5)
ax_pos_err_y.plot(t,-3*sqrt(covr[1,1,:]), 'g--', label='-3 sigma', linewidth = 0.5)
ax_yaw_err.plot(t,yaw-xr[2,:], 'r', label='yaw error', linewidth = 1.0)
ax_yaw_err.plot(t,3*sqrt(covr[2,2,:]), 'g--', label='+3 sigma', linewidth = 0.5)
ax_yaw_err.plot(t,-3*sqrt(covr[2,2,:]), 'g--', label='-3 sigma', linewidth = 0.5)
ax_pos_err_x.legend()
ax_pos_err_y.legend()
ax_yaw_err.legend()

fig_grdtrk_estimate, ax_grdtrk = plt.subplots()
ax_grdtrk.plot(pos[0,:], pos[1,:], 'b', label='ground truth')
ax_grdtrk.plot(xr[0,:], xr[1,:], 'r', label='estimation')
ax_grdtrk.set_xlabel('x')
ax_grdtrk.set_ylabel('y')
ax_grdtrk.legend()

#fig_test_ellipse, ax_te = plt.subplots()
for i in range(0, N, int(N/10)):
    cov_ellipse_from_meancov(xr[0:2, i], covr[0:2, 0:2, i], ax_grdtrk, edgecolor = 'red')

'''
uncertainty propagation without fix (histogram)
'''
Q = np.diag([(T*vel_var)**2, (T*vel_var)**2, (T*w_var)**2])
R = np.diag([lidar_var**2, lidar_var**2])

# set initial condition  
xr = np.zeros((3,N))
covr = np.zeros((3,3,N))
covr[:,:,0] = 1e-10 * np.eye(3)



'''
eskf estimation
'''

'''
hf estimation
'''

plt.show()






