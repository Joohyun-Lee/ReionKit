class CosmoSim:
    def __init__(
        self,
        boxsize=64,
        h=0.7,
        Omega_m=0.3,
        Omega_lambda=0.7,
        Omega_b=0.05,
        sigma8=0.8,
        ns=0.96
    ):
        """
        Initialize the cosmological simulation parameters.
        
        Parameters
        ----------
        - boxsize: Box size [cMpc/h]
        - h: Hubble constant H0/100 [km/s/Mpc]
        - Omega_m: Matter density parameter
        - Omega_lambda: Dark energy density parameter
        - Omega_b: Baryon density parameter
        - sigma8: Normalization of the power spectrum
        - ns: Spectral index of the primordial power spectrum
        """
        self.boxsize = boxsize
        self.h = h
        self.Omega_m = Omega_m
        self.Omega_lambda = Omega_lambda
        self.Omega_b = Omega_b
        self.sigma8 = sigma8
        self.ns = ns

        
    def __repr__(self):
        return (f"Boxsize={self.boxsize}\
                CosmologyParams(h={self.h}, Omega_m={self.Omega_m},\
                Omega_lambda={self.Omega_lambda}, Omega_b={self.Omega_b},\
                sigma8={self.sigma8}, ns={self.ns})")
    