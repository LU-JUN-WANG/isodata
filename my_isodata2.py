import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist, squareform
class my_isodata():
    def __init__(self, data, k0, nmin, dmin, sigma):
        self.__data = np.asarray(data)
        self.__k0 = k0
        self.__nmin = nmin
        self.__dmin = dmin
        self.__sigma = sigma

    def assign(self,d, c):
        dm = cdist(d,c)
        ac = np.argmin(dm, axis=1)
        return ac
    def fit(self):
        #步驟1:隨機挑選中心
        center = self.__data[np.random.choice(len(self.__data), self.__k0, replace=False), :]
        # center = np.insert(center,0,values=range(center.shape[0]),axis=1)
        #步驟2:針對數據集中每個樣本，分配到與某距離最小的中心的類中
        assignation_center = self.assign(self.__data,center[:,1:3])
        # report = np.insert(report,center.shape[1],values=assignation_center,axis=1)
        #步驟3:判斷上述每個類的樣本數是否小於nmin，如果小於，則把該類拋棄並把其元素重新分配到其他類
        while True:
            unique, counts = np.unique(assignation_center, return_counts=True)#找個中心元素數量
            frequency = np.asarray((unique, counts)).T #統計各類元素數量
            if len(np.where(frequency[:,1]<self.__nmin)[0])>0:
                first_nonmin = np.where(frequency[:,1]<self.__nmin)[0][0]#找第一個小於n_min的中心
                center = np.delete(center,first_nonmin,axis=0)
                assignation_center = self.assign(self.__data, center[:, 1:3])
                continue
            else:
                break
        report = np.insert(self.__data, center.shape[1], values=assignation_center, axis=1)
        #步驟4:針對每個類別，重新計算質心
        for i,v in enumerate(np.unique(assignation_center)):
            w = np.where(report[:,-1]==v)
            center[i,:]=np.meaen(self.__data[w,:],axis=0)
        #步驟5:如果當前中心數量少於k0/2則分裂
        if center.shape[0]<self.__k0/2:
            std = np.array([])
            for i in center[0].shape:
                mask = report[:,-1]==i
                choose_point = report[mask,0:2]
                std = np.insert(std,np.std(choose_point,axis=0),axis=0)
            max_std = np.max(std,axis=1)
            max_stddim = np.argmax(std,axis=1)
            where_max = np.where(max_std>self.__sigma)
        #步驟6:如果當前中心數量大於2*k0則合併
        elif center.shape[0]>self.__k0*2:
            while True:
                D = np.triu(squareform(pdist(center)))
                if len(np.where((D==D)&(D!=0)&(D<self.__dmin)))>0:
                    r = np.where((D == D) & (D != 0) & (D <= self.__dmin))
                    first_nodmin = list(map(list, zip(r[0], r[1])))[0]
                    mask = report[:,-1]==first_nonmin[0] | report[:,-1]==first_nodmin[1]
                    tem_var = report[mask,:]
                    report = np.delete(report,mask,axis=0)
                    choose_center = center[first_nodmin,:]
                    center=np.delete(center,first_nodmin,axis=0)
                    six_ac = self.assign(report[:,0:2],center)
                    report = np.insert(report[:,0:2],center.shape[1],values=six_ac,axis=1)
                    new_center = np.average(choose_center,weights=[[tem_var[tem_var[:,-1]==first_nodmin[0]]],[tem_var[tem_var[:,-1]==first_nodmin[1]]]])
                    center = np.insert(center,new_center,axis=0)
                    tem_var2 = tem_var[:,0:2]
                    tem_var2 = np.insert(tem_var2,tem_var2.shape[1],values=center.shape[0],axis=1)
                    report = np.vstack((report,tem_var2))
                else:
                    break
        #步驟7:如果達到最大迭代次數則終止，沒有則回到步驟2










