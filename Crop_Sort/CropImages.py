import glob
import cv2
import matplotlib.pyplot as plt

# Gets all of the images from a folder and stores it in matfiles. 
direc = '/home/jeff/Documents/Ultrascan/Ultrascan_Py/Crop_Sort/Test_Images/'
matfiles = glob.glob(direc + '*')

nfiles = len(matfiles)

# Loops through all the images and makes a cropped version
# that is then written into whereever this file is locoated. 
for i in range(nfiles):
    img_name = matfiles[i].replace(direc, '').replace('.bmp', '')
    img_number = int(img_name.split("_")[1])
    
    img = cv2.imread(matfiles[i])
    crop_img = img[0: 640, 0: 640]
    
    if (img_number < 10):
        cv2.imwrite('image_000' + str(img_number) + '.bmp', crop_img)
    elif (img_number < 100):
        cv2.imwrite('image_00' + str(img_number) + '.bmp', crop_img)
    elif (img_number < 1000):
        cv2.imwrite('image_0' + str(img_number) + '.bmp', crop_img)
    elif (img_number < 10000):
        cv2.imwrite('image_' + str(img_number) + '.bmp', crop_img)