import numpy as np
import skimage.transform
import glob
import cv2
import math
import os

# Variables
crop = 1
adjust = 1

# Instead of hard-coding, implement a file finder
direc = '/home/jeff/Documents/Ultrascan/Ultrascan_Py/Plantar/Scan_3d__bmp/'

padSide = 0
numSlices = 512  # Number of rotational slices
rotTot = 255     # Total rotation sweep in degrees
probeID = 0
imgW = 640       # slice image width
imgH = 480       # slice image height
thic = 3         # Thicknes of slice in pixels
interval = 1     # Interval between slices, 1 => use all slices
imgDepth = 12       # Depth of scan in cm
dia = 16         # Diameter of circle in cm
filtSize1 = 1
filtSize2 = 1

# Opens the file and store the text in a 2d array
direc_to_data_txt = direc + 'Scan_3d__data.txt'

fileID = open(direc_to_data_txt)
data_txt_2d = np.genfromtxt(fileID,dtype='str')
fileID.close()

for i in range(len(data_txt_2d)):
    for j in range(len(data_txt_2d[i])):
        print(data_txt_2d[i][j], end=" ")
    print()
    
# Storing the data read from Scan_3d__data.txt into variables
probeID = int(data_txt_2d[15][1])
imgW = int(data_txt_2d[16][1])
imgH =  int(data_txt_2d[17][1])
imgDepth = int(data_txt_2d[1][1]) / 10
rotTot = int(data_txt_2d[19][1])
dia = float(data_txt_2d[20][1]) * 2

# This is to crop the sides of the image in the conical shape for 7.5 Hz
# b is the value that corresponds to guardband and cone edges
# m is the value that corresponds to the diagonal
b = 70 # these b and m values are the original ones that Austeen had
m = .5
gband = 0

if (probeID == 512):
    if (imgDepth == 4):
        b = 0;
        m = 0
    elif (imgDepth == 8):
        b = 70
        m = 0.5
    elif (imgDepth == 12):  
        b = 0
        m = 0
    elif (imgDepth == 16): 
        b = 0
        m = 0
elif (probeID == 256): 
    if (imgDepth == 2):
        b = 320
    elif (imgDepth == 4):
        b = 200 + gband  
    elif (imgDepth == 6):   
        b = 132.5 + gband; 
    elif (imgDepth == 8):
        b = 100 + gband
    elif (imgDepth == 10): 
        b = 80 + gband

# Taking the first part of the image to crop the top area 
# where it is convex
cropTop = round(imgH / (imgDepth * 4))

# Bottom tangent line to the circle of imaging
# will be used to crop the bottom of the image
cropBottom = imgH - round((imgDepth - dia / 2) * imgH / imgDepth)

# Pixels to pad array for rotation
# so that we do not lose any data while slicing
cropPadSide = round((dia - imgDepth) * imgH / imgDepth)

if (cropPadSide < 0):
    # Only in case of 16cm imaging - need to pad sides to prevent loss
    cropPadSide = abs(cropPadSide)
    padSide = 1
else:
    padSide = 0

# Pixels to take off sides
# Based on the C# code in the software sdk we have settings of 4cm interval
# We can adjust the C# code and get intervals of 2cm as well
imgW = 640 
cropSide_d5 = 0
startPix = 0

'''
if (imgDepth == 4):
    imgW = 640
    startPix = 80 # What is startPix? I don't see it initialized anywhere else
if (imgDepth == 8):
    imgW = 640
    startPix = 200
if (imgDepth == 12):
    imgW = 640
    startPix = 240
if (imgDepth == 16):
    imgW = 640
    startPix = 260
'''

# For our knowledge
if (crop == 1):
    print("Cropping Images: ON")
else:
    print("Cropping Images: OFF")

# For our knowledge
if (adjust == 1):
    print("Adjusting Images: ON")
else:
    print("Adjusting Images: OFF")
    
print("Always adjust the number of silces analyzed and always adjust the rotation angle of the images")

# Places all the file paths to the images in nf
direc_to_bitmaps = direc + 'Scan_3d__bmp/' + 'Image_*_*_*.bmp'
nf = glob.glob(direc_to_bitmaps)

# Sorts nf by image number
image_number = [0] * len(nf)

for i in range(len(nf)):
    image_number[i] = int(nf[i].replace(direc + "Scan_3d__bmp/Image_", "").split("_")[0])
    
image_number = np.array(image_number)
indices = np.argsort(image_number)
nf_placeholder = nf.copy() # Without copy, it would just point back to nf

for i in range(len(indices)):
    nf[i] = nf_placeholder[indices[i]]
    
numSlices = len(nf)

print("Directory: " + direc)
print("Number of slices: " + str(numSlices))
print("Diameter: " + str(dia))
print("Rotation Angle: " + str(rotTot))

rotInc = rotTot / numSlices # Degrees rotated between slices

# Creates a 3d array of zeroes
summed = np.uint8(np.zeros((imgH + cropPadSide, imgW, imgH + cropPadSide)))
[L,W,H]= summed.shape

print("Reconstructing...")

# Creating the slices
for i in range(0, numSlices, interval):
    filename = nf[i]
    rawImg = np.array(cv2.imread(filename, cv2.IMREAD_GRAYSCALE))
    
    # Size of the image before cropping - always 480x640
    SImg1 = rawImg.shape
    height, width = rawImg.shape
    
    # if (adjust == 1):
        # rawImg = imguidedfilter()
        # rawImg = imadjust()
    
    if (crop == 1):
        if (probeID == 512):
            rawImg[0: cropTop - 1, : ] = 0
            rawImg[cropBottom - 1: height - 1, : ] = 0
            # Make a much smaller raw image after cropping
            # the top and the bottom components
            rawImg = rawImg[ : , startPix : startPix + imgW]
    
    # Crop the image into a triangle/cone shape while cutting off the sides
    # Use top width and bottom width and shave off the sides
    # of the image turning it into a polygon
    if (probeID == 512):
        for q in range(0, imgH):
            cropTriangle = round(imgW / 2 - (m * q + b))
            rawImg[q, 0 : cropTriangle - 1] = 0
            rawImg[q, imgW - cropTriangle - 1: imgW - 1] = 0
    elif (probeID == 256):
        for q in range(0, imgH):
            cropSide = round(imgW / 2 - b)
            rawImg[q, 0 : cropSide - 1] = 0
            rawImg[q, imgW - cropSide - 1: imgW - 1] = 0
    
    # Filtering using an averaging method to make brightness - 1
    # for bad areas
    rawImg[0, : ] = -1
    
    # New width and height in case rawImg got cropped
    height, width = rawImg.shape
    
    # Filtering using an averaging method to make brightness -1 for all bad areas
    for k in range(0, width):
        for m in range(0, height - 6):
            if ((rawImg[m + 6 - 1, k - 1] / 7 + rawImg[m + 5 - 1, k - 1] / 7 +
               rawImg[m + 4 - 1, k - 1] / 7 + rawImg[m + 3 - 1, k - 1] / 7 +
               rawImg[m + 2 - 1, k - 1] / 7 + rawImg[m + 1 - 1, k - 1] / 7) > 90):
                break
            rawImg[m - 1, k - 1] = -1
            
    # Size of raw image after the cropping
    SImg2 = rawImg.shape
    
    # Only PadSide for 16cm arrays and greater
    # Initialize 3d array to contain slice
    temp = np.uint8(np.zeros((imgH + cropPadSide, imgW, (2 * thic) + 1)))
    
    for j in range(-thic - 1, thic):
        if (padSide == 0):
            temp[ : , : , thic + 1 + j - 1] = np.vstack((rawImg, np.zeros((cropPadSide, imgW))))
        else:
            temp[ : , : , thic + 1 + j - 1] = np.vstack((np.zeros((cropPadSide, imgW)), rawImg))
            
    # Permute 3d slice, rotate about vertical axis, permute back
    theta = rotInc * (i - 1)
    temp = np.transpose(skimage.transform.rotate(np.transpose(temp, (0, 2, 1)), theta), (0, 2, 1))
    
    L2, W2, H2 = temp.shape
    
    if (L2 > L):
        temp = temp[0: L - 1, : , : ]
        L2 = L
    if (H2 > H):
        temp = temp[ : , : , 0: H - 1]
        
    off1 = math.floor((L - L2) / 2)
    off2 = math.floor((H - H2) / 2)
    
    # Getting the max value from summed[0 + off1: L2 - 1 + off1, : , 0 + off2: H2 - 1 + off2]
    # and placing that value for all summed[0 + off1: L2 - 1 + off1, : , 0 + off2: H2 - 1 + off2]
    max_value = np.amax(summed[0 + off1: L2 - 1 + off1, : , 0 + off2: H2 - 1 + off2])
    summed[0 + off1: L2 - 1 + off1, : , 0 + off2: H2 - 1 + off2] = max_value
    
    # Output slice number
    print("Slice number: " + str(i + 1))
    
# Export the newly created slices
nanInd = np.isnan(summed)
# summed2 = summed.astype('float32') # 32 bit vs 64 bit - less memory used
# summed2[nanInd] = -1000

convert = np.amax(np.amax(np.amax(summed))) / 255

# Creating X,Y, and Z folders and putting the sliced images into those
print("Outputting New Slices...")

# Creating SlicesZ
# This block of code transposes summed before getting the for loops
# Does a 3d tranpose
try:
    os.mkdir(direc + "NewSlicesZ")
except FileExistsError as e:
    print("NewSlicesZ was already created.")
    
summed1 = np.transpose(summed, (0, 2, 1))
L, W, H = summed1.shape

for j in range(0, W):
    filename2 = direc + "NewSlicesZ/" + "Image_" + '{0:03}'.format(j) + ".bmp"
    cv2.imwrite(filename2, summed1[ : , j, : ])

# Creating SlicesX
try:
    os.mkdir(direc + "NewSlicesX")
except FileExistsError as e:
    print("NewSlicesX was already created.")

summed2 = np.transpose(summed, (1, 0, 2))
L, W, H = summed2.shape

for j in range(0, H):
    filename2 = direc + "NewSlicesX/" + "Image_" + '{0:03}'.format(j) + ".bmp"
    cv2.imwrite(filename2, summed2[ : , : , j])

# Creating SlicesY
try:
    os.mkdir(direc + "NewSlicesY")
except FileExistsError as e:
    print("NewSlicesY was already created.")

summed3 = np.transpose(summed, (1, 2, 0))
L, W, H = summed3.shape

for j in range(0, L):
    filename2 = direc + "NewSlicesY/" + "Image_" + '{0:03}'.format(j) + ".bmp"
    cv2.imwrite(filename2, summed3[j , : , :])
    
# Creating SlicesZ
# This block of code transposes summed in the for loop
# Does a 2d transpose
'''
try:
    os.mkdir(direc + "NewSlicesZ")
except FileExistsError as e:
    print("NewSlicesZ was already created.")
    
for j in range(0, W):
    filename2 = direc + "NewSlicesZ/" + "Image_" + '{0:03}'.format(j) + ".bmp"
    cv2.imwrite(filename2, np.transpose(summed[ : , j, : ], (1, 0)))
    
# Creating SlicesX
try:
    os.mkdir(direc + "NewSlicesX")
except FileExistsError as e:
    print("NewSlicesX was already created.")
    
for j in range(0, H):
    filename2 = direc + "NewSlicesX/" + "Image_" + '{0:03}'.format(j) + ".bmp"
    cv2.imwrite(filename2, np.transpose(summed[ : , : , j], (1, 0)))

# Creating SlicesY
try:
    os.mkdir(direc + "NewSlicesY")
except FileExistsError as e:
    print("NewSlicesY was already created.")
    
for j in range(0, L):
    filename2 = direc + "NewSlicesY/" + "Image_" + '{0:03}'.format(j) + ".bmp"
    cv2.imwrite(filename2, np.transpose(summed[j , : , :], (1, 0)))
'''