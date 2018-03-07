from skimage import io
from matplotlib import pyplot as plt
import numpy as np
import cmath
from math import *
def main():
    fileName = "Kolo.jpg"
    #fileName = "CT_ScoutView.jpg"
    #fileName = "xx.jpg"
    image = io.imread("tomograf-zdjecia/" + fileName, as_grey=True)
    h, w = image.shape #wysokosc, szerokosc
    #print(image.shape)
    r = min(h, w) / 2 - 1 # aby sie odsunac od krawedzi
    center = h/2, w/2
    krok = 360
    ndetektorow = 10
    rozpietosc = 30 # phi
    alfa = 360/krok
    #emiter = center[0], center[1] + r
    #print(emiter)
    #image[int(emiter[0])][int(emiter[1])] = 1
    for i in range(0, krok):
        xe = int(round(r * cos(radians(alfa * i)) + center[1],1))
        ye = int(round(r * sin(radians(alfa * i)) + center[0],1))
        image[ye][xe] = 1
        print("emiter=",xe,ye)
        for k in range(0, ndetektorow):
            xdi = int(round(r * cos(radians(alfa + 180 - rozpietosc / 2 + k * rozpietosc / (ndetektorow - 1))), 1) + center[1])
            ydi = int(round(r * sin(radians(alfa + 180 - rozpietosc / 2 + k * rozpietosc / (ndetektorow - 1))), 1) + center[0])
            image[ydi][xdi] = 1
            #print(xdi.real, ydi.real)
    io.imshow(image)
    plt.show()
    #print(image)
    #print(center)
    #print(cos(radians(30)))
if __name__ == '__main__':
    main()