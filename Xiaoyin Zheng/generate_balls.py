import numpy as np
import random
import matplotlib.pyplot as plt
from skimage import io


def ball(img,r):
    '''

    This function is to generate a ball in img with radius r and centers randomly picked.
    '''
    # img_size=z,y,x
    t=img.shape
    len_x=t[2]
    len_y=t[1]
    len_z=t[0]

    centerx=random.randint(int(r),int(len_x-r))
    centery=random.randint(int(r),int(len_y-r))
    centerz=random.randint(int(r),int(len_z-r))
    # print(centerz)
    
    x=np.array(range(len_x)).reshape(1,1,len_x)
    y=np.array(range(len_y)).reshape(1,len_y,1)
    z=np.array(range(len_z)).reshape(len_z,1,1)
    mask=(x-centerx)**2+(y-centery)**2+(z-centerz)**2 <= r**2
    return mask


img=np.zeros((200,200,200))
r=np.random.normal(12,6,100) #list of radius of the ball following gaussian distribution

img0=img.copy()
rl=[]
for i in range(len(r)):
    '''
    Generate 100 balls in img with size 200*300*200
    '''
    mm_t=ball(img0,r[i])
    temp=np.argwhere(mm_t)
    temp=np.array(temp)
    s_overlap=np.argwhere(img[tuple(temp.T)]==1)
    '''
    Ensure there is no overlap for the data
    '''
    while len(s_overlap)>0:
        mm_t=ball(img0,r[i])
        temp=np.argwhere(mm_t)
        temp=np.array(temp)
        s_overlap=np.argwhere(img[tuple(temp.T)]==1)
    
    print(len(s_overlap))
    if len(s_overlap)==0:
        img=img+mm_t
        rl.append(r[i])
        print('{} is ok, number={}'.format(r[i],i))

'696 for data cropping'
big=np.zeros((200,300,200))
big[:,50:250,:]=img



io.imsave('696_100balls_noverlap_12_6.tiff',big)




