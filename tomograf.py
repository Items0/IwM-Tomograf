from skimage import io
from matplotlib import pyplot as plt
import numpy as np
import cmath
from math import *
from skimage.transform import radon

def calculateAvg(myLine, image):
    sum = 0
    for y, x in myLine:
        sum += image[y,x]
    return sum / len(myLine)

def line(x0,y0, x1,y1):
    myLine = []
    dx = x1-x0
    dy = y1-y0

    def sign(x):
        if x >= 0: return +1
        else: return -1

    inc_x = sign(dx) # uwzględnienie znaków dx
    inc_y = sign(dy) # i dy

    dx = abs(dx)     # teraz odcinek został "przeniesiony"
    dy = abs(dy)     # do właściwego oktantu

    if dx >= dy:     # dy/dx <= 1 -- odcinek leży w "niebieskim" oktancie

        d = 2 * dy - dx
        delta_A = 2 * dy
        delta_B = 2 * dy - 2 * dx
        x, y = (0, 0)
        for i in range(dx+1):
            myLine.append([y0 + y, x0 + x])
            if d > 0:
                d += delta_B
                x += inc_x
                y += inc_y
            else:
                d += delta_A
                x += inc_x
    else:            # dy/dx > 1 -- odcinek leży w "czerwonym" oktancie
                     # proszę zwrócić uwagę na wspomnianą zamianę znaczenia zmiennych
        d = 2 * dx - dy
        delta_A = 2 * dx
        delta_B = 2 * dx - 2 * dy
        x, y = (0, 0)
        for i in range(dy+1):
            myLine.append([y0 + y, x0 + x])
            if d > 0:
                d += delta_B
                x += inc_x
                y += inc_y
            else:
                d += delta_A
                y += inc_y
    return myLine

def main():
    fileName = "Kwadraty2.jpg"
    #fileName = "Paski2.jpg"
    #fileName = "CT_ScoutView.jpg"
    #fileName = "Shepp_logan.jpg"
    image = io.imread("tomograf-zdjecia/" + fileName, as_grey=True)
    h, w = image.shape #wysokosc, szerokosc
    r = min(h, w) / 2 - 1 # aby sie odsunac od krawedzi
    center = h / 2, w / 2
    krok = 180
    ndetektorow = 100
    rozpietosc = 270 # phi
    alfa = 360 / krok
    tab = np.zeros((krok, ndetektorow))
    revertTab = np.zeros((h,w))
    #create sinogram
    for i in range(0, krok):
        xe = int(round(r * cos(radians(alfa * i)) + center[1], 1))
        ye = int(round(r * sin(radians(alfa * i)) + center[0], 1))
        for k in range(0, ndetektorow):
            xdi = int(round(r * cos(radians(alfa * i + 180 - rozpietosc / 2 + k * rozpietosc / (ndetektorow - 1))), 1) + center[1])
            ydi = int(round(r * sin(radians(alfa * i + 180 - rozpietosc / 2 + k * rozpietosc / (ndetektorow - 1))), 1) + center[0])
            tab[i][k] = calculateAvg(line(xe, ye, xdi, ydi), image)

    #restore image from sinogram
    for i in range(0, krok):
        xe = int(round(r * cos(radians(alfa * i)) + center[1], 1))
        ye = int(round(r * sin(radians(alfa * i)) + center[0], 1))
        for k in range(0, ndetektorow):
            xdi = int(round(r * cos(radians(alfa * i + 180 - rozpietosc / 2 + k * rozpietosc / (ndetektorow - 1))), 1) + center[1])
            ydi = int(round(r * sin(radians(alfa * i + 180 - rozpietosc / 2 + k * rozpietosc / (ndetektorow - 1))), 1) + center[0])
            myLine = line(xe, ye, xdi, ydi)
            for y,x in myLine:
                revertTab[y][x] += tab[i][k]

    #normalization
    mymax = np.max(revertTab)
    for i in range(0, revertTab.shape[0]):
        for k in range(0, revertTab.shape[1]):
            revertTab[i][k] = revertTab[i][k] / mymax

    io.imshow(revertTab)
    plt.show()

if __name__ == '__main__':
    main()