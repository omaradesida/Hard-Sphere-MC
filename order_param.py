def alt_nematic_order_param(configs,Nbeads,Nchains):
    
    """
    
    Plots the nematic order parameter over time for an array of Atoms Objects or a trajectory.
    
    configs: Array of Atoms Objects/Trajectory for the atoms objects which needs to be calculated
    Nbeads: Number of beads in the chain which is being calculated
    Nchains: The number of chains present in the system which is being calculated
    
    """

    order_param_array=[]
    Nconfigs = len(configs) 


    for img in configs:
        
        Q_sum = np.zeros((3,3)) # Initialise Q_sum
        
        for i in range(Nchains):
            
            # Calculate the vector from the beginning chain to the end chain

            a0 = i*Nbeads      # index of zeroth bead on chain i
            a1 = i*Nbeads + 3  # index of 4th bead on chain i 

            # Find end-to-end vector
            chainvector = img.get_distance(a0, a1, mic=True, vector=True)
            
            
            # Normalise chainvector
            chainvector = chainvector/np.linalg.norm(chainvector)

            # Accumulate Q_sum
            Q_sum = Q_sum + 1.5 * np.outer(chainvector, chainvector) - 0.5 * np.identity(3)
            
        print(Q_sum)
        # Find the average Q
        Qav = Q_sum/Nchains
        print(Qav)
        # Diagonalise this
        eigenvals, eigenvects = np.linalg.eig(Qav)
        
        print("eigen")
        print(eigenvals)
        print(eigenvects)
        
        # The maximum positive eigenvalue is the P2 order parameter
        P2 = np.amax(eigenvals)
        
        order_param_array.append(P2)

        
        # If we wanted the director then this is the corresponding eigenvector
        #k = np.where(eigenvals == P2)
        #director = eigenvects[k]
        
        
    
    return order_param_array
