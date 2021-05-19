#Author Name: Varun Kankanallu
#Class: CME696 Spring

from pystackreg import StackReg
from skimage import io, transform
from matplotlib import pyplot as plt
from scipy import stats
import numpy as np

#load the two images you want to align from your folder location
path = "C:/Users/K R Varun/Python_project/"
one = path + "BottomCenter_Ti_resemble_0001.tiff"
two = path + "BottomCenter_Ti_resemble_0002.tiff"

image1 = io.imread(one)
image2 = io.imread(two)

#view the two images you want to align

fig, axes = plt.subplots(1, 2, figsize=(8, 4))
ax = axes.ravel()
ax[0].set_xlabel('Grid_values_along_x')
ax[0].set_ylabel('Grid_values_along_y')

ax[1].set_xlabel('Grid_values_along_x')
ax[1].set_ylabel('Grid_values_along_y')

ax[0].set_title('120x120 Grid view of image1')
ax[1].set_title('120x120 Grid view of image2')

ax[0].imshow(image1)
ax[1].imshow(image2)

plt.savefig('20210518_Grid_view_image1_and_image2.png', format='png', dpi=300)
fig.tight_layout()
plt.show()

#print the array of each image, this is so that we can better understand the image we are trying to align
print(image1)
print(image2)

#to understand the characteristics of the array we are working with
print(np.amax(image1))
print(image1.ndim)
print(image1.shape)
print(np.amax(image2))

#normalised pixel intensity plot along the x and y-axis
#left to right
cols1 = image1.mean(axis=0)
ncols1 = cols1/np.amax(cols1)
maxpc1 = np.where(cols1 == np.amax(cols1))
print(maxpc1[0])
# bottom to top
rows1 = image1.mean(axis=1)
nrows1 = rows1/np.amax(rows1)
maxpr1 = np.where(rows1 == np.amax(rows1))

# plot histograms
f, ax = plt.subplots(2, 1)
ax[0].plot(ncols1)
ax[1].plot(nrows1)
f.show()

# left to right
cols2 = image2.mean(axis=0)
ncols2=cols2/np.amax(cols2) #Normalisation of the array to the maximum

# bottom to top
rows2 = image2.mean(axis=1)
nrows2=rows2/np.amax(rows2) #Normalisation of the array to the maximum

#offset calculation evaluated by the
maxpr2 = np.where(rows2 == np.amax(rows2)) #Pointer to the maximum value of the rows array
maxpc2 = np.where(cols2 == np.amax(cols2)) #Pointer to the maximum value of the columns array
print(maxpc1[0]-maxpc2[0]) #to obtain the x-coordinates of the Offset
print(maxpr1[0]-maxpr2[0]) #to obtain the y-coordinates of the Offset

# plotting the two plots together so that we can visually see the drift between the two images
f, ax = plt.subplots(1,2, figsize=(14,4))
ax[0].plot(ncols2, label="image1")
ax[1].plot(nrows2, label="image1")
ax[0].plot(ncols1, label="image2")
ax[1].plot(nrows1, label="image2")

ax[0].legend(loc='lower right', shadow=True, fancybox=True)
ax[1].legend(loc='upper right', shadow=True, fancybox=True)

ax[0].set_title('Normalised pixel intensity along x-axis for image1 and image2')
ax[1].set_title('Normalised pixel intensity along Y-axis for image1 and image2')

plt.savefig('20210518_Normalised_PI_along_x_y_for_image1_and_image2.png', format='png', dpi=300)

f.show()

fig = plt.figure()
ax = fig.add_subplot(111)

L=stats.linregress(ncols1[0:115],ncols2[0:115])
print(L)

ax.plot(ncols1[0:115], L.slope*ncols1[0:115]+L.intercept, color='green', linestyle='dashed', marker='o', markerfacecolor='green', markersize=2)

ax.scatter(ncols1,ncols2, label= "Image_1 vs. Image_2, Along x-axis")
ax.legend(loc='lower right', shadow=True, fancybox=True)
ax.set_xlabel('Image1 values')
ax.set_ylabel('Image2 values')
ax.set_title('linear correlation between Image1 and Image2')
plt.show()

from skimage.feature import register_translation
shifted, erro, diffphase=register_translation(image1, image2, 100)
print(f"Detected subpixel offset (y,x): {shifted}")

from scipy.ndimage import shift
corrected_image=shift(image2, shift=(1.19,-4.18), mode='constant')

fig, axes = plt.subplots(1, 2, figsize=(8, 4))
ax = axes.ravel()
ax[0].set_xlabel('Grid_values_along_x')
ax[0].set_ylabel('Grid_values_along_y')

ax[1].set_xlabel('Grid_values_along_x')
ax[1].set_ylabel('Grid_values_along_y')

ax[0].set_title('120x120 Grid view of image1')
ax[1].set_title('120x120 Grid view of Corrected_image')

ax[0].imshow(image1)
ax[1].imshow(corrected_image)

plt.savefig('20210518_Grid_view_image1_and_corrected_image.png', format='png', dpi=300)
fig.tight_layout()
plt.show()

#correlation plot between corrected_image and image_1
#left to right
cols3 = corrected_image.mean(axis=0)
ncols3=cols3/np.amax(cols3) #Normalisation of the array to the maximum
maxpc3 = np.where(cols3 == np.amax(cols3)) #Pointer to the maximum value of the array
print(maxpc3[0])

# bottom to top
rows3 = corrected_image.mean(axis=1)
nrows3=rows3/np.amax(rows3)
maxpr3 = np.where(rows3 == np.amax(rows3))

# plot histograms
f, ax = plt.subplots(1,2, figsize=(14,4))
ax[0].plot(ncols3)
ax[1].plot(nrows3)

ax[0].set_title('Normalised pixel intensity along x-axis for Corrected_image')
ax[1].set_title('Normalised pixel intensity along y-axis for Corrected_image')

plt.savefig('20210518_Normalised_PI_along_x_y_for_Corrected_image.png', format='png', dpi=300)
f.show()

fig = plt.figure()
ax = fig.add_subplot(111)

L1=stats.linregress(ncols1[0:115],ncols3[0:115])
print(L1)
ax.plot(ncols1[0:115], L.slope*ncols1[0:115]+L.intercept, color='green', linestyle='dashed', marker='o',
     markerfacecolor='green', markersize=2)
ax.scatter(ncols1[0:115],ncols3[0:115], label= "Corrected_image vs. Image_1, Along x-axis")
ax.legend(loc='lower right', shadow=True, fancybox=True)
ax.set_ylabel('Corrected_image values')
ax.set_xlabel('Image_1 values')
ax.set_title('linear correlation between Corrected_image and Image_1')
plt.savefig('20210518_Correlation_plot_Corrected_image_vs_Image_1.png', format='png', dpi=300)
plt.show()


#Correlation plot between offset_image and image_1
from scipy.ndimage import shift
offset_image=shift(image2, shift=(1,-5), mode='constant')

#left to right
cols4 = offset_image.mean(axis=0)
ncols4=cols4/np.amax(cols4) #Normalisation of the array to the maximum
maxpc3 = np.where(cols4 == np.amax(cols4)) #Pointer to the maximum value of the array
print(maxpc3[0])

# bottom to top
rows4 = offset_image.mean(axis=1)
nrows4=rows4/np.amax(rows4)
maxpr4 = np.where(rows4 == np.amax(rows4))

fig = plt.figure()
ax = fig.add_subplot(111)

L1=stats.linregress(ncols1[0:115],ncols4[0:115])
print(L1)
ax.plot(ncols1[0:115], L.slope*ncols1[0:115]+L.intercept, color='green', linestyle='dashed', marker='o',
     markerfacecolor='green', markersize=2)
ax.scatter(ncols1[0:115],ncols3[0:115], label= "Offset vs. Image_1, Along x-axis")
ax.legend(loc='lower right', shadow=True, fancybox=True)
ax.set_ylabel('Corrected_image values')
ax.set_xlabel('Image_1 values')
ax.set_title('linear correlation between offset_image and Image_1')
plt.savefig('20210518_Correlation_plot_offset_image_vs_Image_1.png', format='png', dpi=300)
plt.show()




