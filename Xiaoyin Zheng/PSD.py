from scipy import ndimage
import numpy as np
import matplotlib.pyplot as plt
from skimage import io
import time


img=io.imread('696_100balls_noverlap_12_6.tiff')
'''
Crop the data according to the peak found previously
'''
img_crop=img[:,63:240,:]
dis_filled2 = ndimage.morphology.distance_transform_edt(img_crop)
dis_filled_crop=dis_filled2.copy()
print(dis_filled_crop.shape)
# Find all the pixels have particle size more than 20 and assign a value of -1


#200*300*200
# xlen = 200
# ylen = 300
# zlen = 200

x0=0
x1=200
y0=63
y1=240
z=200


rs_l=[]                 # The list of radius
volume_culmulative=[]

for i in range(1,30,2):
    rs_l.append(i)

for rs in rs_l:
    t1=time.time()
    dis_filled_centers=dis_filled_crop.copy()

    '''
    Select the pixels with their value larger than rs in the distance map
    '''
    for i in range(0,z):
        for j in range(0,y1-y0):
            for k in range(0,x1-x0):
                if dis_filled_centers[i,j,k]>=rs:
                    dis_filled_centers[i,j,k]=-1
    
    '''
    Assign these pixels with a value 1
    '''
    dis_filled_centers_inv=dis_filled_centers.copy()
    for i in range(0,z):
        for j in range(0,y1-y0):
            for k in range(0,x1-x0):
                if dis_filled_centers_inv[i,j,k] != -1:
                    dis_filled_centers_inv[i,j,k]=0
                else:
                    dis_filled_centers_inv[i,j,k]=1
                    
    '''Mask'''
    mask11=np.zeros((z,y1-y0,x1-x0))
    
    mask11= img_crop==1
    
    
    struct2 = ndimage.generate_binary_structure(3, 1)
    iteration=70
    dis_filled_centers_inv=ndimage.binary_erosion(dis_filled_centers_inv, structure=struct2, iterations=1).astype(dis_filled2.dtype)
    '''
    dilation and calculate the total volume with size larger than rs
    '''
    d1=ndimage.binary_dilation(dis_filled_centers_inv,structure=struct2,mask=mask11,iterations=iteration).astype(dis_filled2.dtype)
    vt=sum(sum(sum(d1)))
    volume_culmulative.append(vt)
    print(vt)
    t2=time.time()
    print('{} is finished. Time used:{}'.format(rs,t2-t1))


'''
Plot the particle size distribution
'''
volume=[]
for i in range(len(rs_l)-1):
    volume.append(volume_culmulative[i]-volume_culmulative[i+1])

volume=np.array(volume)
volume=volume/sum(volume)
x=rs_l[:-1]
plt.scatter(x,volume)
plt.plot(x,volume)
plt.xlabel('size/pixel')
plt.ylabel('volume fraction')
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
