# -*- coding: utf-8 -*-
"""
@author: Vahid Haseltalab
"""
from stl import mesh
import numpy as np
from scipy import spatial
import shapely.geometry as sg


### Slicing an STL file format. It requires STL address or vertices of a stl file,
### the height of slice plane and direction of the slicing (xy,xz,yz)
### Slicing an STL file format. It requires STL address, and the height of slice plane
class Slice:
    def vertices(self,a):
        # importing ths stl file
        msh = mesh.Mesh.from_file('%s'%a)
        zvalues = [] ### stores the z value of a vertex 
        ### categorizing the vertices into a list based on their faces
        vrt = [[] for s in range(len(msh))] 
        for i in range(len(msh)):
            p1 = (msh[i][:3]).tolist()
            p2 = (msh[i][3:6]).tolist()
            p3 = (msh[i][6:]).tolist()
            zvalues.append(p1[-1]);zvalues.append(p2[-1])
            zvalues.append(p3[-1])
            vrt[i] = [p1,p2,p3]
        return vrt,zvalues
    
    def __init__(self,addr = 'bunny.stl',direction = 'xy'):
        if type(addr) == str:
            self.vrt,self.zvalues = self.vertices(addr)
            self.direction = direction
        else:
            self.direction = direction
            self.vrt = addr
            if direction == 'xy':
                self.zvalues = addr.reshape([len(addr)*3,3])[:,2]
            elif direction == 'xz':
                self.zvalues = addr.reshape([len(addr)*3,3])[:,1]
            elif direction == 'yz':
                self.zvalues = addr.reshape([len(addr)*3,3])[:,0]
     ### finds intersection coordinates between a point and a line 
    def eqn1(self,p1,p2,z):
        if p1[2]==p2[2]: return tuple([p1[0],p1[1]])
        ### the ratio that can apply to distance of x and y 
        t = (z-p1[2]) / (p2[2]-p1[2]) 
        return (p1[0] + (p2[0]-p1[0])*t , p1[1] + (p2[1]-p1[1])*t)
    def eqn2(self,p1,p2,z):
        if p1[1]==p2[1]: return tuple([p1[0],p1[2]])
        ### the ratio that can apply to distance of x and y 
        t = (z-p1[1]) / (p2[1]-p1[1]) 
        return (p1[0] + (p2[0]-p1[0])*t , p1[2] + (p2[2]-p1[2])*t)
    def eqn3(self,p1,p2,z):
        if p1[0]==p2[0]: return tuple([p1[1],p1[2]])
        ### the ratio that can apply to distance of x and y 
        t = (z-p1[0]) / (p2[0]-p1[0]) 
        return (p1[1] + (p2[1]-p1[1])*t , p1[2] + (p2[2]-p1[2])*t)
    ### checks whether the z plane is crossing through the line
    def checkline(self,zl,z):
        if z <= np.max(zl) and z >= np.min(zl):
            return True
        else: return False
    ### finds intersection coordinates between a plane and a triangular facet        
    def trintersct(self,l,z):
        l = np.array(l)
        if self.direction == 'xy' and (l[:,2] == z).all(): return []
        elif self.direction == 'xz' and (l[:,1] == z).all(): return []
        elif self.direction == 'yz' and (l[:,0] == z).all(): return []
        inlst = []
        for i in range(3):
            pt1 = l[i];pt2 = l[i-1]
            if self.direction == 'xy':
                zl = [pt1[2],pt2[2]]
                if self.checkline(zl,z):
                    p = self.eqn1(pt1,pt2,z)
                    inlst.append(p)
            elif self.direction == 'xz':
                zl = [pt1[1],pt2[1]]
                if self.checkline(zl,z):
                    p = self.eqn2(pt1,pt2,z)
                    inlst.append(p)
            elif self.direction == 'yz':
                zl = [pt1[0],pt2[0]]
                if self.checkline(zl,z):
                    p = self.eqn3(pt1,pt2,z)
                    inlst.append(p)
        if len(inlst) == 3:
            if np.array_equal(np.round(inlst[0],4),np.round(inlst[1],4)): del(inlst[1])
            else: del(inlst[2])
        if inlst:
            if np.array_equal(np.round(inlst[0],4),np.round(inlst[1],4)): return []
        return inlst
    #The function clear will remove empty lists from S
    def clear(self,A):
        n=len(A);s=0;k=0
        for i in range(n):
            s+=int(not A[i])
        W=[[0 for i in range(1)]for j in range(n-s)]
        for i in range(n):
            if A[i]:
                W[k]=A[i]
                k+=1
        return W
    def totuple(self,a): ## convert an array to nested tuple
        try:
            return tuple(self.totuple(i) for i in a)
        except TypeError:
            return a
    def order(self,pnts,prec = 5):
        if not pnts: return []
        pnts = np.round(pnts,prec)
        pnts = self.totuple(pnts)
        wires = [];Q1 = set()
        for pn in pnts:
            #pn = simplify(pn,0.01)
            if pn not in Q1:
                Q = list(pn)
                Q1.add(pn)
                a1 = Q[0]; a2 = Q[1]; a3 = Q[-2]; a4 = Q[-1]
                c = 0
                while c < 1:
                    run = False
                    for b in pnts:
                        #b = simplify(b,0.01)
                        if b not in Q1:
                            b1 = b[0]; b2 = b[-1]
                            if a4 == b1:
                                run = True
                                Q = Q + list(b)
                                Q1.add(b)
                                a3 = b1; a4 = b2
                            elif a1 == b2:
                                run = True
                                Q = list(b) + Q
                                Q1.add(b)
                                a1 = b1; a2 = b2
                            elif a4 == b2:
                                run = True
                                Q = Q + list(b)[::-1]
                                Q1.add(b)
                                a3 = b2; a4 = b1
                            elif a1 == b1:
                                run = True
                                Q = list(b)[::-1] + Q
                                Q1.add(b)
                                a1 = b2; a2 = b1
                    if not run: c+=1
                wires.append(Q)
        if len(wires) == 1:
            return wires[0]
        if prec == 0: return wires[0]
        else: return self.order(wires,prec-1)
    def cross_section(self,z):
        z1 = np.array(self.zvalues).reshape((int(len(self.zvalues)/3),3))
        idx = np.where((z1>z-5)&(z1<z+5))[0]
        fnum = set()
        L = []
        for a in list(set(idx)):
            if a not in fnum:
                fnum.add(a)
                cr = self.trintersct(self.vrt[a],z)
                if cr: L.append(cr)
        L1 = np.round(L,3)
        lines = list(so.polygonize(L1.tolist()))
        if lines:
            return list(lines[0].exterior.coords)
        else:
            return list(self.order(L))
