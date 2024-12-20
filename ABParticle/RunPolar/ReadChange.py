import subprocess
def FindKeyPosition(f,keyword,p1):
    f.seek(p1)
    data = f.read(len(keyword))
    while data != keyword:
        data = f.read(len(keyword))
        p1 = f.tell()
        f.seek(p1+1 - len(keyword))
    return p1

def ReadChange(targetdir,keywords,UseWords):
    with open(targetdir,"rt+") as f:
        p1 = f.tell()
    # print(p1)
        for i in range(len(keywords)):
            p1 = FindKeyPosition(f,keywords[i],p1)
            f.seek(p1)
            f.write(UseWords[i])
            f.flush()
    return 0
            
  
  
  
basedir = "./"
datadir = basedir
names = ["../src/field/ABParticle/ABParticleFieldClass.cu","bacteriaPolar.cu","create_matrix.py","DrawTurbulence.py","visual_video.py","Bash.sh"]
# AA = -20
# NumParticles = 800
for name in names:
    # print(name)
    # keywords = ["AA","NumParticles"]
    # usewords = [AA,NumParticles]
    ReadChange(datadir+name,["test"],["06"])
# names1 = []
# for name in names1:
#     print(name)
#     ReadChange(datadir+name,[""])
subprocess.run(['bash',basedir+"Bash.sh" ])

            
            
        
        
              
# example
   
# basedir = "/home/lyuan/share/LLC/T_phase/"

# a_array = ["2.000","1.900","1.800","1.700","1.600",
#            "1.500","1.400","1.300","1.200","1.100",
#            "1.000","0.900","0.800","0.700","0.600",
#            "0.500","0.400","0.300","0.200","0.100",
#            "0.000","-0.10","-0.20","-0.30","-0.40",
#            "-0.50","-0.60","-0.70","-0.80","-0.90",
#            "-1.00","-1.10","-1.20","-1.30","-1.40"]

# c_array = ["0.01","0.02","0.03","0.04","0.05","0.06",
#            "0.07","0.08","0.09","0.10","0.11","0.12",
#            "0.13","0.14","0.15","0.16","0.17","0.18",
#            "0.19","0.20","0.21","0.22","0.23"
#            ]
# for a in a_array:
#     for c in c_array:
#         keywords = ["23/"," a =","cmin1 ="]
#         UseWords = ["a"+a+"c"+c,a,c]
#         datadir = basedir + "data/"+"CoolPhase1023/"+UseWords[1]
#         if os.path.exists(datadir) ==0:
#             os.makedirs(datadir)
            
            
#         ReadChange(datadir+"bacteria_test.cu",keywords,UseWords)
        
        
        
        
        
# with open(basedir+"bacteria_test.cu","rt+") as f:
#     # for i in range(line1):
#     #     data = f.readline()
#     p1 = f.tell()
#     # print(p1)
#     for i in range(3):
#         p1 = FindKeyPosition(f,keywords[i],p1)
#         f.seek(p1)
#         f.write(UseWords[i])
#         f.flush()
        
        
        
        
        
        
        
        
        
    # print(data)
    # p1 = 3116
    # f.seek(p1)
    # data = f.read(3)
    
    # while (data != "23/"):
    #     data = f.read(3)
    #     p1 = f.tell()
    #     # print(f"p1={p1}")
    #     # print(f"data = {data}")
    #     f.seek(p1-2)
    # f.seek(p1)
    # # f.write("a"+a+"c"+c)
    # f.flush()
    # while (data != " a ="):
    #     data = f.read(4)
    #     p1  = f.tell()
    #     # print(p1)
    #     f.seek(p1-3)
    # f.seek(p1)
    # f.write(a)
    # f.flush()
    # # print(f"{a}")
    
    # while (data!= "cmin1 ="):
    #     data = f.read(7)
    #     # print(p1)
    #     p1 = f.tell()
    #     f.seek(p1-6)
    # f.seek(p1)
    # f.write(c)
  
        
    
 
    
