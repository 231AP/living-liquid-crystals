import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import cv2
import os  
# savenames = [
#     "videoa-2c0", "videoa-2c001", "videoa-2c002", "videoa-2c003", "videoa-2c004", "videoa-2c005",  "videoa-2c01" ,"videoa-2c02", "videoa-2c03", "videoa-2c04", "videoa-2c05", "videoa-2c1", "videoa-2c2" ,"videoa-2c3", "videoa-2c4", "videoa-2c5",
#     "videoa-05c0", "videoa-05c001", "videoa-05c002", "videoa-05c003", "videoa-05c004", "videoa-05c005",  "videoa-05c01" ,"videoa-05c02", "videoa-05c03", "videoa-05c04", "videoa-05c05", "videoa-05c1", "videoa-05c2" ,"videoa-05c3", "videoa-05c4", "videoa-05c5",
#     "videoa0c0", "videoa0c001", "videoa0c002", "videoa0c003", "videoa0c004", "videoa0c005",  "videoa0c01" ,"videoa0c02", "videoa0c03", "videoa0c04", "videoa0c05", "videoa0c1", "videoa0c2" ,"videoa0c3", "videoa0c4", "videoa0c5",
# "videoa-1c0", "videoa-1c001", "videoa-1c002", "videoa-1c003", "videoa-1c004", "videoa-1c005",  "videoa-1c01" ,"videoa-1c02", "videoa-1c03", "videoa-1c04", "videoa-1c05", "videoa-1c1", "videoa-1c2" ,"videoa-1c3", "videoa-1c4", "videoa-1c5",
#  "videoa1c0","videoa1c001", "videoa1c002", "videoa1c003", "videoa1c004", "videoa1c005",  "videoa1c01" ,"videoa1c02", "videoa1c03", "videoa1c04", "videoa1c05", "videoa1c1", "videoa1c2" ,"videoa1c3", "videoa1c4", "videoa1c5",
# "videoa05c0", "videoa05c001", "videoa05c002", "videoa05c003", "videoa05c004", "videoa05c005",  "videoa05c01" ,"videoa05c02", "videoa05c03", "videoa05c04", "videoa05c05", "videoa05c1", "videoa05c2" ,"videoa05c3", "videoa05c4", "videoa05c5",
# "videoa2c0", "videoa2c001", "videoa2c002", "videoa2c003", "videoa2c004", "videoa2c005",  "videoa2c01" ,"videoa2c02", "videoa2c03", "videoa2c04", "videoa2c05", "videoa2c1", "videoa2c2" ,"videoa2c3", "videoa2c4", "videoa2c5",]
savenames = [
    "test11"
]
for savename in savenames:
    datadir ="../photo_video/" + savename +"/"
    savedir ="../photo_video/" + savename + "/"
    videoname = savedir + "../" + savename + ".mp4"
    image_folder = savedir 

    image_files = [f for f in os.listdir(image_folder) if f.endswith('.png')]  

 
    image_files.sort()  
    print(image_files)
    frame = cv2.imread(os.path.join(image_folder, image_files[0]))
    height, width, layers = frame.shape
 
    fourcc = cv2.VideoWriter_fourcc(*'mp4v')  
    

    out = cv2.VideoWriter(videoname, fourcc, 10.0, (height, width))  
    

    for image_file in image_files:  
        img = cv2.imread(os.path.join(image_folder, image_file))  
        out.write(img)  
    

    out.release()
    