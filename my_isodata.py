import numpy as np
from scipy.spatial.distance import pdist, squareform
class my_isodata():
    def __init__(self, data, k0, nmin, dmin, sigma):
        self.__data = data
        self.__k0 = k0
        self.__nmin = nmin
        self.__dmin = dmin
        self.__sigma = sigma
    def dest(self,pdata,center_list):
        center = []
        for i in center_list:
            center.append(locals()[i])
        min_index = np.argmin(squareform(pdist(pdata, center)), axis=1)
        for i,c in enumerate(min_index):
            locals()[center_list[c]].append(pdata[i])
        for i in center_list:
            return locals()[i]
    def id_merge(self,center_list):
        center = []
        for i in center_list:
            center.append(locals()[i])
        d = np.triu(pdist(center,center))
        a = np.where((d~=0)&(d<self.__dmin))

    def fit(self):
        fcenter = np.random.choice(self.__data,self.__k0,replace=False)
        fmin_index = np.argmin(squareform(pdist(self.__data,fcenter)),axis=1)
        totalc=[]
        for i in range(self.__k0):
            locals()['c'+str(i)]=[]
            totalc.append('c'+str(i))
        tem_var = []
        for i in range(self.__k0):
            if locals()['c'+str(i)]<self.__nmin:
                tem_var=.append(locals()['c'+str(i)])
                del locals()['c'+str(i)]


