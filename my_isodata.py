import numpy as np
from scipy.spatial.distance import cdist, squareform
class my_isodata():
    def __init__(self, data, k0, nmin, dmin, sigma):
        self.__data = data
        self.__k0 = k0
        self.__nmin = nmin
        self.__dmin = dmin
        self.__sigma = sigma
    def fit(self):
        center = self.__data[np.random.choice(len(self.__data),self.__k0, replace=False), :]
        dt = cdist(self.__data,center)
        cofdic = dict(zip(range(len(center)),center))
        totalc = list(range(len(center)))
        belong = np.nanargmin(dt,axis=1)
        for i in np.unique(belong):
            locals()['c'+str(i)] = []
        for i in range(len(self.__data)):
            locals()['c'+str(belong[i])].append(self.__data[i])
        del_c = []
        for c in totalc:
            if len(locals()['c'+str(c)]) < self.__nmin:
                del_c.append(c)
                center[c] = None
                td = cdist(locals()['c'+str(c)], center)
                tb = np.nanargmin(td,axis=1)
                for w in range(len(tb)):
                    tw = tb[w]
                    locals()['c' + str(tw)].append(locals()['c' + str(c)][w])
        ntotalc = list(set(totalc)-set(del_c))
        for c in ntotalc:
            m = np.mean(locals()['c'+str(c)],axis=0)
            center[c] = m
        if len(ntotalc)<self.__k0/2:
            d = np.triu(squareform(pdist(center)))
            wd = np.where((d==d)&(d!=0)&(d<self.__dmin))
