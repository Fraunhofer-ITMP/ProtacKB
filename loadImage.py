#from Ptac_V2 import createGraph

#a=createGraph()
#print(len(a))

from PIL import Image
import glob
image_list = []
path = r"C:\Users\reagon.karki\PycharmProjects\ProtacKB\data\test_images\*.png"
for file in glob.glob(path):
    im=Image.open(file)
    image_list.append(im)


print(len(image_list))
print(image_list[0])