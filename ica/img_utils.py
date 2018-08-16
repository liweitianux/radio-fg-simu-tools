import numpy
import scipy
import mdp
import sys
import math
import exceptions


def is_peak(img,i,j):
    if i<1 or j<1 or i>=img.shape[0]-1 or j>=img.shape[1]-1:
        return False

    if img[i+1,j]>=img[i,j]:
        return False
    
    if img[i-1,j]>=img[i,j]:
        return False

    if img[i,j+1]>=img[i,j]:
        return False

    if img[i,j-1]>=img[i,j]:
        return False

    if img[i+1,j+1]>=img[i,j]:
        return False

    if img[i-1,j+1]>=img[i,j]:
        return False
    
    if img[i+1,j-1]>=img[i,j]:
        return False

    if img[i-1,j-1]>=img[i,j]:
        return False

    return True

def is_bottom(img,i,j):
    if i<1 or j<1 or i>=img.shape[0]-1 or j>=img.shape[0]-1:
        return False

    if img[i+1,j]<=img[i,j]:
        return False
    
    if img[i-1,j]<=img[i,j]:
        return False

    if img[i,j+1]<=img[i,j]:
        return False

    if img[i,j-1]<=img[i,j]:
        return False

    if img[i+1,j+1]<=img[i,j]:
        return False

    if img[i-1,j+1]<=img[i,j]:
        return False
    
    if img[i+1,j-1]<=img[i,j]:
        return False

    if img[i-1,j-1]<=img[i,j]:
        return False

    return True


def make_round_region(img,center_i,center_j,r):
    result=[]
    for i in range(int(center_i-r),int(center_i+r+2)):
        for j in range(int(center_j-r),int(center_j+r+2)):
            if i>=0 and i<img.shape[0] and j>=0 and j<img.shape[1]:
                if (i-center_i)**2+(j-center_j)**2<=r**2:
                    result.append((i,j))

    return result

def make_box_region(img,center_i,center_j,h,w):
    result=[]
    for i in range(int(center_i-h/2),int(center_i+h/2)):
        for j in range(int(center_j-w/2),int(center_j+w/2)):
            result.append((i,j))

    return result


def make_pie_region(img,center_i,center_j,ri,ro,a1,a2):
    result=set()
    
    center_x=center_j
    center_y=center_i
    r=ri
    while r<=ro:
        angle=a1
        while angle<=a2:
            ii=center_i+r*math.sin(angle/180.*math.pi)
            jj=center_j+r*math.cos(angle/180.*math.pi)
            if ii>=0 and jj>=0 and ii<img.shape[0] and jj<img.shape[1]:
                
                result.add((ii,jj))
            angle+=.1/r/math.pi*180

        r+=.5
    return list(result)




def enum_peaks(img):
    result=[]
    for i in range(0,img.shape[0]):
        for j in range(0,img.shape[1]):
            if is_peak(img,i,j):
                result.append((i,j))

    return result

def max_element(img):
    result=None
    max_value=-1E99
    for i in range(0,img.shape[0]):
        for j in range(0,img.shape[1]):
            if max_value<img[i,j]:
                max_value=img[i,j]
                result=(i,j)
    return result

def draw_circle(ci,cj,r):
    i=r
    j=0
    F=0
    n=r
    result=set([])

    while n>0:
        result.add((ci+i,cj+j))
        result.add((ci+j,cj+i))
        result.add((ci-i,cj-j))
        result.add((ci-j,cj-i))
        result.add((ci-i,cj+j))
        result.add((ci-j,cj+i))
        result.add((ci+i,cj-j))
        result.add((ci+j,cj-i))
        if F<=0:
            F+=2*j+1
            j+=1
        else:
            F+=1-2*i
            i-=1
            n-=1
        if i==j:
            result.add((ci+i,cj+j))
    return list(result)
