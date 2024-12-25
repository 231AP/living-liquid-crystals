from moviepy.editor import ImageSequenceClip
from PIL import Image
import os
def png_to_mp4(png_folder, output_file, fps):

    images = [os.path.join(png_folder, img) for img in sorted(os.listdir(png_folder)) if img.endswith(".png")]
    
    if not images:
        raise ValueError("no PNG")
    first_image = Image.open(images[0])
    width, height = first_image.size  
    clip = ImageSequenceClip(images, fps=fps)
    clip = clip.resize(newsize=(width, height))
    clip.write_videofile(output_file, codec='libx264')
savenames = [
    "test34","test33"
    # "new1"
]
for savename in savenames:
    png_folder = "../photo_video/"+ savename + "/"
    output_file = png_folder + "../" + savename + ".mp4"
    # output_file = png_folder + "../" + "new1" + ".mp4"

    fps =10  




    png_to_mp4(png_folder, output_file, fps)
