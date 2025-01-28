from moviepy import ImageSequenceClip
from PIL import Image
import os
def png_to_mp4(png_folder, savedir, fps):
    output_file = savedir + ".mp4"
    images = [os.path.join(png_folder, img) for img in sorted(os.listdir(png_folder)) if img.endswith(".png")]
    
    if not images:
        raise ValueError("no PNG")
    first_image = Image.open(images[0])
    width, height = first_image.size  
    clip = ImageSequenceClip(images, fps=fps)
    # clip = clip.resize(newsize=(width, height))
    clip.write_videofile(output_file, codec='libx264')

if __name__ == "__main__":
    savenames = [
        "test93",
        # "new1"
    ]
    for savename in savenames:
        png_folder = "../photo_video/"+ savename + "/"
        output_file = png_folder + "../" + savename 
        # output_file = png_folder + "../" + "new1" + ".mp4"

        fps =10  




        png_to_mp4(png_folder, output_file, fps)
