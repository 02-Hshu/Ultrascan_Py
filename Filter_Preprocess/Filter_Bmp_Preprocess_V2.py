import numpy as np
import cv2
import matplotlib.pyplot as plt
import glob
import os

# Removes the blackspots from the img.
def remove_blackspots(img):
    # Tranposes Img1 then flattens by rows.
    # Same thing as flattening Img1 by columns.
    Img1_V = np.transpose(img, (1, 0)).flatten()
    
    # Returns indices of the zero elements.
    # Turns all the 0s to 128s, essentially removes the black spots.
    Zeros_Img1 = np.where(np.array(Img1_V) == 0)[0]
    
    for i in range(0, len(Zeros_Img1)):
        j = Zeros_Img1[i]
        Img1_V[j] = 128
    
    [length, width] = img.shape
    
    Img1_V_Gray = np.reshape(Img1_V, (width, length))
    Img1_V_Gray_true = np.transpose(Img1_V_Gray, (1, 0))
    
    return Img1_V_Gray_true

# g(x) = a * f(x) + b
# a and b are gain and bias parameters, said to control
# constrast and brightness respectively. 

# Increases the contrast of the img.
def increase_contrast(img, alpha):
    # np.clip clips all the values to fit between the set interval
    # Using np.array is faster than using 2 nested for loops to
    # change each element. 
    img = np.clip(np.array(img) * alpha, 0, 255)
    
    return img

# Increases the sharpness of the img.
def increase_sharpness(img):
    # Kernel for sharpening
    kernel = np.array([[ 0, -1,  0],
                       [-1,  5, -1],
                       [ 0, -1,  0]])
    
    # Sharpen the image
    img = cv2.filter2D(img, -1, kernel)
    
    return img

# Gets all the bmps and places it in a variable.
direc = '/home/jeff/Documents/Ultrascan/Ultrascan_Py/Plantar/Scan_3d__bmp/'
nf = glob.glob(direc + 'Scan_3d__bmp/Image_*_*_*.bmp')

# Creates a folder to store the newly processed images.
try:
    os.mkdir(direc + 'Preprocessed_Images')
except FileExistsError as e:
    print("Folder was already created")

# Path to newly created folder. 
path_to_preprocessed_images = direc + 'Preprocessed_Images/'

# Loops through all the bmps, applies the image processing functions,
# then writes the newly processed images. 
for i in range(len(nf)):
    # Print out the slice number
    print("Slice number: " + str(i + 1))
    
    # Remove the blackspots from the img
    img = cv2.imread(nf[i], cv2.IMREAD_GRAYSCALE)
    img = remove_blackspots(img)   
    
    # Increase sharpness
    img = increase_sharpness(img)
    
    # Increase constrast
    img = increase_contrast(img, 1.4)
    
    # Writes to the folder containing the folder with all the scans
    filename = nf[i].replace(direc, '').replace('Scan_3d__bmp/', '')
    cv2.imwrite(path_to_preprocessed_images + filename, img)