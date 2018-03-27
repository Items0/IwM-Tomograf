from skimage import io
from matplotlib import pyplot as plt
import numpy as np
from math import *


def normalizeM(nowyTab):
    nowyTabMax = np.max(nowyTab)
    for i in range(0, nowyTab.shape[0]):
        for k in range(0, nowyTab.shape[1]):
            nowyTab[i][k]=nowyTab[i][k] / nowyTabMax
            if(nowyTab[i][k]<0):
                nowyTab[i][k]=0
    return nowyTab

def normalize(tab):
    mymax = np.max(tab)
    mymin=np.min(tab)
    newMax=1
    newMin=0
    for i in range(0, tab.shape[0]):
        for k in range(0, tab.shape[1]):
            tab[i][k] = (((tab[i][k]-mymin)/(mymax-mymin))*(newMax-newMin))+newMin
    print(np.min(tab))
    print(np.max(tab))
    return tab

def createMask(ndetektorow):
    ndetektorow = int(ndetektorow//18) #wielkosc maski = dlugosc wiersza / 4
    if(ndetektorow%2==0):
        ndetektorow+=1
    maska=np.zeros((ndetektorow))
    iterator=0
    for h in range(int(len(maska)/2),len(maska),1):
        if(iterator % 2 != 0):
            maska[h] = -4 / (pi * pi * iterator*iterator)
        else:
            maska[h]=0
        iterator+=1
    print(maska)
    for h in range(0,int(len(maska)/2)):
        y=((-1*h)+len(maska))-1
        maska[h]=maska[y]
    srodek = int(len(maska)/2)
    maska[srodek]=1
    print(maska)
    return maska

'''
def myConvolve(krok, ndetektorow, tab, maska):
    wyn = np.zeros((krok, ndetektorow + len(maska)))
    nowyTab = np.zeros((krok, ndetektorow + len(maska))) #rozszerzona

    for i in range(krok): #skopiowana
        for k in range (len(maska), len(maska) + ndetektorow):
            nowyTab[i][k] = tab[i][k-len(maska)]

    for i in range(krok):
        for k in range(ndetektorow + len(maska)):
            for el in range(len(maska)):
                wyn[i][k] += nowyTab[i][k-el] * maska[el]
    wyn = wyn[:, len(maska):]
    return wyn
'''
def bladSredniokwadratowy(pocz, wynikowa):
    print("pocz",pocz.shape[0],pocz.shape[1])
    print("wynikowa",wynikowa.shape[0],wynikowa.shape[1])
    pocz=pocz[wynikowa.shape[0]:wynikowa.shape[1]]
    return np.power(np.subtract(wynikowa, pocz), 2).mean() ** 0.5

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

def myConvolve2(tab,mask):
    out_sinogram = np.zeros_like(tab)
    base_sinogram=tab
    filtr = mask

    half_filtr_len = len(filtr) // 2
    padding = 2 * half_filtr_len + 1
    pad_array = np.zeros(half_filtr_len)
    for ind_x, x in enumerate(base_sinogram):
        for i in range(base_sinogram.shape[1]):
            out_sinogram[ind_x][i] = np.sum(np.concatenate([pad_array, x, pad_array])[i:i + padding] * filtr)
    return out_sinogram


def main():
    #createMask(301)
    #fileName = "Kolo.jpg"
    fileName = "Kwadraty2.jpg"
    #fileName = "Paski2.jpg"
    #fileName = "CT_ScoutView.jpg"
    #fileName = "Shepp_logan.jpg"
    #fileName="SADDLE_PE.JPG"
    image = io.imread("tomograf-zdjecia/" + fileName, as_grey=True)
    h, w = image.shape #wysokosc, szerokosc
    r = min(h, w) / 2 - 1 # aby sie odsunac od krawedzi
    center = h / 2, w / 2
    krok = 360
    ndetektorow = 301
    rozpietosc = 250 # phi
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

    mask = createMask(ndetektorow)

    #splot

    #afterConv = myConvolve(krok, ndetektorow, tab, mask)
    #afterConv=normalizeM(afterConv)
    #afterConv = np.zeros((krok, ndetektorow))
    #afterConv=[]
    '''
    tmp=[]
    for a in range(krok):
        tmp=np.convolve(tab[a],mask,mode='same')
        for i in range(len(tmp)):
            afterConv[a][i]=tmp[i]
        tmp=[]
    print(np.min(afterConv))
    afterConv=normalize(afterConv)
    '''

    afterConv = myConvolve2(tab,mask)

    for i in range(0, krok):
        xe = int(round(r * cos(radians(alfa * i)) + center[1], 1))
        ye = int(round(r * sin(radians(alfa * i)) + center[0], 1))
        for k in range(0, ndetektorow):
            xdi = int(round(r * cos(radians(alfa * i + 180 - rozpietosc / 2 + k * rozpietosc / (ndetektorow - 1))), 1) + center[1])
            ydi = int(round(r * sin(radians(alfa * i + 180 - rozpietosc / 2 + k * rozpietosc / (ndetektorow - 1))), 1) + center[0])
            myLine = line(xe, ye, xdi, ydi)
            for y, x in myLine:
                revertTab[y][x] += afterConv[i][k]
    revertTab = normalize(revertTab)
    io.imshow(tab)
    plt.show()
    io.imshow(afterConv)
    plt.show()
    io.imshow(revertTab)
    plt.show()

if __name__ == '__main__':
    main()