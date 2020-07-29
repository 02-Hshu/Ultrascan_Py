import glob
import os
import cv2
import matplotlib.pyplot as plt

# Gets all of the images from a folder and stores it in matfiles. 
# Make sure to change these when selecting a new folder with images to crop 
direc = '/home/jeff/Documents/Ultrascan/Training_Images/TS2_Data/'
direc_to_images = direc + '/Images/'
direc_to_cropped_images = direc + '/Cropped_Images/'

# Create cropped images folder
try:
    os.mkdir(direc_to_cropped_images)
except OSError:
    print("Failed to create cropped images directory or cropped images directory already created")

# Make sure this leads to diretory with images
matfiles = glob.glob(direc_to_images + '*')

nfiles = len(matfiles)

# Loops through all the images and makes a cropped version
# that is then written into whereever this file is locoated. 
for i in range(nfiles):
    # Make sure it replaces the path to the directory with the images
    img_name = matfiles[i].replace(direc_to_images, '').replace('.bmp', '')
    img_number = int(img_name.split("_")[1])
    
    img = cv2.imread(matfiles[i])
    crop_img = img[0: 640, 0: 800]
    
    if (img_number < 10):
        cv2.imwrite(direc_to_cropped_images + 'image_000' + str(img_number) + '.png', crop_img)
    elif (img_number < 100):
        cv2.imwrite(direc_to_cropped_images + 'image_00' + str(img_number) + '.png', crop_img)
    elif (img_number < 1000):
        cv2.imwrite(direc_to_cropped_images + 'image_0' + str(img_number) + '.png', crop_img)
    elif (img_number < 10000):
        cv2.imwrite(direc_to_cropped_images + 'image_' + str(img_number) + '.png', crop_img)
