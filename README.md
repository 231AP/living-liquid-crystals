# living-liquid-crystals
date: Jan 28 2025 

除夕快乐！

Hello everyone! I'm Li. Yuan, and Here is the readme code of living-liquid-crystals.

Today I upload a new file named ABP-LLC, it will be the new main file of this project.

Now let me introduct this program and tell you how to use it.

## Introduction

This program is a simulation code of a biophysical system named living-liquid-crystal, actually we are trying to simulate how self-propelled particles move in anisotropy hydrodymic system. 
So we have two parts in our code, one of them is the PDE-Solver, to calculate the hydrodymic background, for example, the velocity field and the direction of liquid crystals, this PDE-Solver's main part was created by my supervisor Prof.You. And the other part of the code is used to be a MD simulation which was created by Benchang Wu and Qizheng Lai, both of them are undergraduate students in Xiamen University. I tryed to fix them to create this new code which could simulation the reaction between fields and particles, using a simple coarse-grained way to make our particles system became a field that can reaction with other filed in PDE-Solver. At the same time, the particles system could read the data from field system and tell every particle how to move. 

## How to use

Here is the simplest way to use my code to do some simulation and data analysis.
First of all, making sure you are working on a computer with Nvidia's GPU, and you need to install CUDA on your computer. Actually, we suggest that you could just use cloud sever to do this.

After that, python is needed, as in my code, I'd like to use python to contral and draw photos and videos. And you need to try to have this environment to make sure you can use this could:

numpy(pip install numpy), pandas(pip install pandas), moviepy(pip install moviepy), cv2 (pip install opencv-python), matplotlib(pip install matplotlib).

After you have python, you can use the commands in brackets to install the environment we need.

And now, you can open the file you download from this project, open the LLC-ABP and get into the file named RunPolar, try use this command:

python main.py

if all the preparation is OK, the code will start and all the data export will be saved in the file /ABP-LLC/data, photos and videos will be saved in the file /ABP-LLC/photo_video.

After that, if you want to change parameters into whatever you want, you can change it in the main.py, and you can find all the paremeters I haved created in the class LLCparams, in the LLC-ABP.py.

while you have calculated too many simulation, you can go to the analysis file, there are some codes may help you to do some data analysis.

After all, if you have any problem and advise, it will be welcome to email 770395058@qq.com, actually add the qq 770395058 will be more quickly.

	
