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
        ### the ratio that can apply to distance of x and y 
        t = (z-p1[2]) / (p2[2]-p1[2]) 
        return [p1[0] + (p2[0]-p1[0])*t , p1[1] + (p2[1]-p1[1])*t]
    def eqn2(self,p1,p2,z):
        ### the ratio that can apply to distance of x and y 
        t = (z-p1[1]) / (p2[1]-p1[1]) 
        return [p1[0] + (p2[0]-p1[0])*t , p1[2] + (p2[2]-p1[2])*t] 
    def eqn3(self,p1,p2,z):
        ### the ratio that can apply to distance of x and y 
        t = (z-p1[0]) / (p2[0]-p1[0]) 
        return [p1[1] + (p2[1]-p1[1])*t , p1[2] + (p2[2]-p1[2])*t] 
    ### checks whether the z plane is crossing through the line
    def checkline(self,zl,z):
        if z < max(zl) and z > min(zl):
            return True
        else: return False
    ### finds intersection coordinates between a plane and a triangular facet        
    def trintersct(self,l,z):
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
    
    def optlst(self,L2,p,r):
        idx1 = np.where((L2[:,0]>p[0]-r)&(L2[:,0]<p[0]+r))[0]
        L3 = L2[idx1]
        L4 = L3[L3[:,1].argsort()]
        idx2 = np.where((L4[:,1]>p[1]-r)&(L4[:,1]<p[1]+r))[0]
        L5 = L4[idx2]
        return L5
    ### Ordering the points in a correct sequence
    def order1(self,L,r):
        l = [[] for _ in range(5)]
        i = 0
        p = L[0]
        del(L[0])
        T1 = spatial.cKDTree(L)
        l[i].append(p)
        while L:
            ds, idx = T1.query(p)
            pp = L[idx]
            if ds > r:
                i += 1
            l[i].append(pp)
            del(L[idx])
            if not L:
                break
            T1 = spatial.cKDTree(L)
            p = pp
        l = self.clear(l)
        for a in l:
            a.append(a[0])
        return l
    def order2(self,L):
        L = np.round(L,5); lst = []
        lst.extend(L[:1].tolist()); a = np.array(lst[-1])
        L = np.delete(L,0,axis=0)
        while L.any():
            vec = L-a; nn = np.linalg.norm(vec,axis=1)
            idx1 = nn.argsort()[:1][0]; lst.append(L[idx1].tolist())
            L = np.delete(L,idx1,axis=0); a = np.array(lst[-1])
        lst.append(lst[0])
        return lst
    ### order the sequence of points to form polygon
    def order3(self,L):
        convex = sg.MultiPoint(L).convex_hull.buffer(3).buffer(-3)
        line = sg.LineString(L)
        k1 = [L.copy() for b in range(3)]
        try: ff = sg.Polygon(L); pl0 = L
        except: pl0 = []
        if convex.length*2 < line.length:
            pl1 = list(convex.exterior.coords)
        else: pl1 = []
        try:
            pl2 = self.order2(k1[1])
            if len(pl2) == 1 and len(pl2[0] > 2): pl2 = pl2[0] 
        except: pl2 = []
        area = [sg.Polygon(a).area for a in [pl0,pl1,pl2]]
        idx = area.index(max(area))
        pl = [pl0,pl1,pl2][idx]
        if len(pl) > 2:
            p = sg.Polygon(pl).buffer(1).buffer(-1)
            if p.area != 0:
                return list(p.exterior.coords)
            else: return list(convex.exterior.coords)
        else: return []
    def cross_section(self,z):
        z1 = np.array(self.zvalues).reshape((int(len(self.zvalues)/3),3))
        idx = np.where((z1>z-4)&(z1<z+4))[0]
        fnum = set()
        L = []
        for a in idx:
            if a not in fnum:
                fnum.add(a)
                L.extend(self.trintersct(self.vrt[a],z))
        pl = self.order3(L)
        if list(pl):
            return pl
        else: 
            print('No intersection!')
            return []