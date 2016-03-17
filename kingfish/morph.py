# An object to store and operate on galaxy morphologies

from numpy import pi, exp, log10, log, gamma

#mu_d = mu_0 + 1.086*(r/Re)

#mu_b = mu_e + b_n *( (r/Re)**(1/n) -1 )
#I_b = I_o np.exp( -(r/Ralpha)**(1/n) )
#L_b_tot = 2*np.pi*Ralpha^2*n*np.gamma(2*n)   
#M_b = mu_0
#I_b = I_e *np.exp( -b_n*((R/Rhalf)**(1/n) -1) )

class Morph(object):

    bulge = {'Ro':None,    #scale length or break radius
             'Re': None,   #half light radius
             'Rhalf':None, # half light radius
             'n':None,
             'b_n':None,
             'So':None,   #mu at scale length
             'Shalf':None,
             'So': None,  
             'Ie':None,
             'Io':None,
             'Ihalf':None,
             'Ltot':None,
             'Lhalf':None,
             'Mhalf':None,
             'Mtot':None} 
    
    disk = {'Re':None,    #scale length
            'Rhalf':None, #half-light radius
            'Se':None,    #mu at scale length
            'Shalf':None, #mu at half -light
            'So':0.       #mu at center
            'Ie':None,
            'Io':None,
            'Ihalf':None,
            'Ltot':None,
            'Lhalf':None,
            'Mhalf':None,
            'Mtot':None} 
    
    
    def __init__(self,  n = None, Shalf = None Ro = None, Re = None, Rhalf = None,):
        pass




class MorphObserved(Morph):
    def _fill_all_parameters(self):
        
        pass


    def mag_within(self, radius = None, Component = 'bulge'):
        pass

    def mu_at(self, radius = None, Component = 'bulge'):
        pass
    
class MorphPhysical(Morph):
    def _fill_all_parameters(self):
        self.bulge['Shalf'] = self.bulge['Se']
        n = self.bulge['n']
        self.bulge['b_n'] = (2.0*n - 1/3. + (4./405.)/n +
                             (46.0/25515.0)/(n**2) + (131.0d/1148175.0)/(n**3))
        self.bulge['So'] = self.bulge.['Se']-1.086*self.bulge['b_n']
        self.bulge['Ro'] = self.bulge['Re']/(self.bulge['b_n']**self.bulge['n'])


    def L_within(self, radius = None, Component = 'bulge'):
        pass

    def I_at(self, radius = None, Component = 'bulge'):

