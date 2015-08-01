import numpy as np

ions = ['NI','NII',
        'OI', 'OII', 'OIII',
        'NeIII', 'NeIV', 'NeV',
        'SII', 'SIII',
        'ClII', 'ClIII', 'ClIV',
        'ArIII', 'ArIV', 'ArV']


class Atom(object):

    cc = np.zeros([5,5])
    qq = np.zeros([5,5])
    jl = np.zeros([5,5])

    def line_ratio(self, tem, den):
        #load temperature insensitive atomic data if necessary
        #load temperature sensitive atomic data if the temperature has changed
        if tem is not self.tem
            self.tem = tem
            self.partii(ion, tem)
            self.hbeta = self.hb_emissivity(tem)
        #solve for line ratio
        return self.emissivity(tem, den)/self.hbeta

    def solve(self, tem, den):
        kt = (1.3807e-16)*tem
        hc = 1.9865e-08
        cq = 8.629e-06
        ts = np.sqrt(tem)

        dat = self.dat
        
        # Transition energy differnces
        e = ( hc * (dat['t'].reshape([5,1])-dat['t'].reshape([1,5])) ).clip(0)
        # Collision Transition probabilities
        q1 = (cq/ts) * (dat['Omega']/dat['weight'].reshape([1,5])) * np.exp(-e/kt)
        q2 = (cq/ts) * (dat['Omega'].T/dat['weight'].reshape([1,5]))
        q = np.trix(q1.T)+np.trix(q2)
        # Critical density 
        ncrit = dat['A'].sum(cumulative =True, axis = 1)
        # matrix manipulation
        a1 = den * q [+a]

        #finally, level populations, emissivities

    def hb_emissivity(self):
        pass

    def emissivity(self):
        pass
        
    def __init__(self):
        self.dat = None
        self.tem = -99
        self.dat = {'A':np.zeros([5,5]),  #radiative transition probability matrix
                    'E':np.zeros(5),      #energy of the level (relative to ground)
                    'Omega':np.zeros([5,5]),
                    'w':np.zeros(5)} # statistical weight of level


class NI(Atom):

    def atomic_data(self):

        a[1:,0] = np.array([7.27e-6, 2.02e-5, 2.71e-3, 6.58e-3])
        a[2:,1] = np.array([1.27e-8, 3.45e-2, 6.14e-2])
        a[3:,2] = np.array([5.29e-2, 2.76e-2])
        e[1:] = np.array([1.92245e-4, 1.92332e-4, 2.88389e-4,2.88393e-4])
        w = np.array([4, 6, 4, 2, 4])

        self.dat['A'] = a
        self.dat['E'] = e
        self.dat['w'] = w
        
    def partii(self, ion, tem):
        t4 = 1.0e-4*tem
        c = np.zeros([5,5])


class OIII(Atom):

    def atomic_data(self):

        a = np.array([[0., 2.62e-5, 3.02e-11, 2.74e-6, 0.],
                     [0., 0., 9.76e-5, 6.74e-3, 2.23e-1],
                     [0., 0., 0., 1.96e-2, 7.85e-4],
                     [0., 0., 0., 0., 1.78],
                     [0., 0., 0., 0., 0.]])
        
        e = np.array([0., 1.132e-6, 3.062e-6, 2.02733e-4, 4.318577e-4])
        w = np.array([1, 3, 5, 5, 1])
        
        self.dat['A'] = a.T
        self.dat['E'] = e
        self.dat['w'] = w
        
    def collision(self, ion, tem):
        t4 = 1.0e-4*tem
        c = np.zeros([5,5])

        cj = np.array([1.835+t4*(0.3981-t4*0.06), 0.2127+t4*(0.0767-0.013*t4)])
        c[1,0] = 0.4825+t4*(0.0806-t4*0.022)
        c[2,0] = 0.2397+t4*(0.0381-t4*0.007)
        c[3:4,0] = cj/9
        c[2,1] = 1.1325+t4*(0.203-t4*0.05)
        c[3:4,1] = cj/3
        c[3:4,2] = cj*5/9
        c[4,3] = 0.3763+t4*(0.3375-t4*0.105)
        
        self.dat[ion]['Omega'] = c

    def getPopulations(self, tem, den, product=True):
        """
        Return array of populations at given temperature and density.

        The method returns a 1-, 2- or 3-D array containing the
        population of each level  for all temperatures and densities
        specified in the input vectors tem and den (which can be
        n-element or 1-element vectors).
        
        If either quantity (tem or den) is a 1-element vector -that
        is, a single value-,  the resulting population array is
        collapsed along that dimension;  as a result, the result
        population array can be a 1-D, 2-D or 3-D array  (the three
        cases corresponding to situations in which both tem and den
        are single values;  one of them is a single value and the
        other an n-element vector; or both are multielement  vectors,
        respectively). In the general case, the level index is the
        first  [WARNING! It is not in physical unit, i.e. ground level
        = 0; to be normalized],  followed by the temperature index (if
        it exists) and the density index.

        Usage:
            O3.getPopulations(1e4, 1e2)
            tem=np.array([10000., 12000., 15000., 20000]) # An array of four temperatures
            den=np.array([600., 800., 1000])      # An array of three densities
            O3.getPopulations(tem, den)           # is a (6, 4, 3) array
            O3.getPopulations(tem, den)[0,2,1]    # Returns the population of level 1 for T = 15000 
                                                    and Ne = 800
            tem = 20000                           # tem is no longer an array
            O3.getPopulations(tem, den)[0,2,1]  # Crashes: one index too much
            O3.getPopulations(tem, den)[0,1]    # Returns the population of level 1 for T = 20000 
                                                    and Ne = 800 [see warning]
            tem=np.array([10000., 15000., 20000]) # An array of three temperatures
            O3.getPopulations(tem, den, product = False)# is a (6, 3) array, tem and den beeing 
                                                            taken 2 by 2.
        
        Parameters:
            - tem       electronic temperature in K
            - den       electronic density in cm^-3
            - product   operate on all possible combinations of temperature and density 
                      (product = True, default case) or on those resulting from combining 
                      the i-th value of tem with the i-th value of den (product = False).
                      If product = False, then tem and den must be the same size.

        """
        tem = np.asarray(tem)
        den = np.asarray(den)
        if product:
        else:
            if tem.shape != den.shape:
                self.log_.error('tem and den must have the same shape', calling=self.calling)
                return None
            res_shape1 = [self.NLevels]
            res_shape_rav1 = [self.NLevels, tem.size]
            res_shape_rav2 = [self.NLevels, self.NLevels, tem.size]
            for sh in tem.shape:
                res_shape1.append(sh)
            tem_rav = tem.ravel()
            den_rav = den.ravel()
            n_level = self.NLevels
            q = self.getCollRates(tem_rav)
            A = self._A
            pop_result = np.zeros(res_shape_rav1)
            coeff_matrix = np.ones(res_shape_rav2)
            sum_q_up = np.zeros(res_shape_rav1)
            sum_q_down = np.zeros(res_shape_rav1)
            sum_A = A.sum(axis=1)
            n_tem = tem_rav.size
            # Following line changed 29/11/2012. It made the code crash when atom_nlevels diff coll_nlevels
            #Atem = np.outer(self._A, np.ones(n_tem)).reshape(n_level, n_level, n_tem)
            Atem = np.outer(self._A[:n_level, :n_level], np.ones(n_tem)).reshape(n_level, n_level, n_tem)
            self._critDensity = Atem.sum(axis=1) / q.sum(axis=1)

            for i in range(1, n_level):
                for j in range(i + 1, n_level):
                    sum_q_up[i] = sum_q_up[i] + q[i, j]
                for j in range(0, i):
                    sum_q_down[i] = sum_q_down[i] + q[i, j]
            for row in range(1, n_level):
                # upper right half            
                for col in range(row + 1, n_level):
                    coeff_matrix[row, col] = den_rav * q[col, row] + A[col, row]
                # lower left half
                for col in range(0, row):
                    coeff_matrix[row, col] = den_rav * q[col, row]
                # diagonal
                coeff_matrix[row, row] = -(den_rav * (sum_q_up[row] + sum_q_down[row]) + sum_A[row])

            vect = np.zeros(n_level)
            vect[0] = 1.
            
            for i in range(tem.size):
                try:
                    pop_result[:, i] = np.linalg.solve(np.squeeze(coeff_matrix[:, :, i]), vect)
                except np.linalg.LinAlgError:
                    pop_result[:, i] = np.nan
                except:
                    self.log_.error('Error solving population matrix', calling=self.calling)
            return np.squeeze(pop_result.reshape(res_shape1))
