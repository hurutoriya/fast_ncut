import numpy as np


class FastNcut:
    def __init__(self, A =  np.array([[2, 0, 0],[0, 1, 0],[0, 0, 0]]) , const=np.array([[0,1]]), max_iter= 1000, opt_tol=1e-12):
        self.A = A
        self.const = const
        self.max_iter = max_iter
        self.opt_tol = opt_tol
        self.B, self.c = self._init_const(self.const)

    def _init_const(self, const):
        const_num = len(const)
        B = np.zeros([const_num,len(self.A)])
        for i, pair in enumerate(const):
            B[i][pair] = [1, -1]
            print(f'{pair[0]+1}th element  and {pair[1]+1}th element be same value in the eigen vector.')
        c = np.zeros([const_num,1])
        return B, c

    def _init_P(self):
        return np.eye(len(self.A))- self.B.T@np.linalg.pinv(self.B@self.B.T)@self.B

    def _init_n(self):
        return self.B.T@ np.linalg.pinv(self.B@self.B.T)@self.c

    def _init_ganma(self,n):
        return np.sqrt(1- np.linalg.norm(n)**2)


    def _init_v(self,ganma, P, n):
        return ganma* P@self.A@n/np.linalg.norm(P@self.A@n)+n


    def _projected_powermethod(self):
        """zeros
        Projected Powermethod optimize a problem
        max_v v^T A v, ||v|| = 1, Bv=c

        Simple constrain: you want to same value i-th and j-th element in eigen vector.
        You can insert 1 and -1 to i-th and j-th element.
        """

        P = self._init_P()
        k = 0
        n_0 = self._init_n()
        ganma = self._init_ganma(n_0)
        # if c is zero vector, escape NaN. 
        if np.count_nonzero(self.B) == 0:
            v = self._init_v(ganma, P, n_0)
        else :
            v = np.random.rand(len(self.A))[np.newaxis].T
        obj = v.T@self.A@v
        obj_old = obj
        print('***Projected Power Method***')

        while  k < self.max_iter:
            v = v/np.linalg.norm(v)
            u = ganma * P@self.A@v/np.linalg.norm(P*self.A*v)
            v = u+n_0
            k+=1
            obj = v.T@self.A@v
            print(f' Iteration : {k} residual : {abs(obj-obj_old)[0][0]}')
            if self.opt_tol > abs(obj-obj_old):
                break
            obj_old = obj
        return v, k

    def fit(self,X):
        const_num = 1
        const_eig_vec, iter_num = self._projected_powermethod()
        print(f'constrained eigen vector : \n {const_eig_vec}')