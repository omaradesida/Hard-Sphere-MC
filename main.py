from joblib import Parallel, delayed
import joblib
import os


def alkane_MC(pressure, solid, maindir, Nbeads = 4, Nchains = 250):
    
    
    
    import numpy as np  # Numpy
    import hs_alkane.alkane as mdl   # Fortran library we'll used to model the alkane chain
    import os 
    

    if solid:
        state='solid'
    else:
        state='fluid'

    workingdir = f"{maindir}{state}P{pressure}"

    if not os.path.exists(workingdir):
        os.makedirs(workingdir) 


    

    #Defining function for creating ase atoms object
    
    from ase import Atoms
    from ase import io

    def mk_ase_config(ibox, Nbeads, Nchains):
        'Uses the current state of the alkane model to construct an ASE atoms object'
    
        # Create and populate ASE object
        model_positions = np.empty([Nchains*Nbeads, 3])
        cell_vectors = mdl.box_get_cell(ibox)
   
        for ichain in range(0, Nchains):
            model_positions[Nbeads*ichain:Nbeads*ichain+Nbeads] = mdl.alkane_get_chain(ichain+1, ibox)
    
        confstring = "C"+str(Nbeads*Nchains)
    
        box_config = Atoms(confstring, positions=model_positions*(1.5/0.4), pbc=True, cell=cell_vectors*(1.5/0.4)) 

        # If we want to wrap all beads inside the box
        box_config.wrap()
    
        return box_config  # Returns ASE atom object



    # Initialise the simulation box and alkane module ready to hold chains

    if os.path.exists(f"{workingdir}/log{state}P{pressure}.log"):
        os.remove(f"{workingdir}/log{state}P{pressure}.log")


    logfile = open(f"{workingdir}/log{state}P{pressure}.log", "w")

    logfile.write(f"Attempting to initialise simulation box")
    print(f"ATTEMPTING TO INITIALISE BOX AT {pressure}")


    mdl.box_initialise()             
    mdl.alkane_set_nbeads(Nbeads)    
    mdl.alkane_set_nchains(Nchains)  
    mdl.alkane_initialise()

    print(f"SUCCESSFULLY INITIALISED BOX AT {pressure}")

    
    logfile.write("system initialised \n")
 
    if os.path.exists(f"{workingdir}/equil.txt"):
        os.remove(f"{workingdir}/equil.txt")


    equilfile = open(f"{workingdir}/equil.txt", "w")  

    



    print(f"GENERATING STRUCTURE AT {pressure}")



    if solid:
        mdl.io_read_xmol() 
    else:


        boxsize = 30.0  # Initial size of cubic box
        ibox = 1
        cell_vectors = np.array([[boxsize,0.0,0.0],[0.0,boxsize,0.0],[0.0,0.0,boxsize]])
        mdl.box_set_cell(ibox, cell_vectors)

        print(f"Setting box size at {pressure}")

        # Don't use link cells when making new chains
        mdl.box_set_bypass_link_cells(1)
        mdl.box_set_use_verlet_list(0)

        print(f"box size set at {pressure}")
        growth_traj = io.Trajectory(f"{workingdir}/growingchains.traj", mode='a')

        ncopy = Nchains
        for ichain in range(ncopy):
            mdl.alkane_set_nchains(ichain+1)  # Only calculate using this and previous chains
            rb_factor = 0.0
            while rb_factor < 0.001:
                print(f"attempting to grow trial chain {ichain+1} at {pressure}")
                rb_factor, ifail = mdl.alkane_grow_chain(ichain+1, ibox,1)
                print(f"trial chain {ichain+1} grown at{pressure}")

                ####
                #growth_traj.write(atoms = mk_ase_config(ibox,Nbeads,Nbeads))
                ####


                if ifail!=0:
                    print(ifail)
                    rb_factor = 0.0

        print(f"RELINKING CELLS AT {pressure}")

        #link cells
        mdl.box_set_link_cell_length(1.6)
        mdl.box_set_bypass_link_cells(0)
        mdl.box_construct_link_cells(ibox)
        mdl.alkane_construct_linked_lists(ibox)
    print(f"STRUCTURE HAS BEEN GENERATED AT {pressure}")

    ibox = 1
    for ichain in range(0, Nchains):
    
    # Checks if internal geometry of chain (bonds/angles) are consistent with the model
        geom_flag = mdl.alkane_check_chain_geometry(ichain+1, ibox)
        if geom_flag != 0:
            logfile.write(f"Incorrect geometry for chain {ichain} for pressure {pressure} in state {state} \n")
            return
       

    # Checks if beads on any two chains overlap
    overlap_flag = mdl.alkane_check_chain_overlap(ibox)
    if overlap_flag != 0:
        logfile.write(f"Overlaps found between chains in initial configuration for pressure {pressure} in state {state} \n")
        return
    else:
        logfile.write(f"No overlaps between chains found in initial configuration for pressure {pressure} in state {state}\n")

    ##Setting up input/output for traj files

    #io.write(f"{workingdir}/configurations.dcd", mk_ase_config(ibox, Nbeads, Nchains), format="dcd")
    #io.write(f"{workingdir}/configurations.dcd", mk_ase_config(ibox, Nbeads, Nchains), format="psf")
    traj = io.Trajectory(f"{workingdir}/configs.traj", mode='a')
    traj.write(atoms = mk_ase_config(ibox,Nbeads,Nchains))
    

    


    if solid:
        mdl.alkane_set_dr_max(0.012)
        mdl.alkane_set_dt_max(0.07)
        mdl.alkane_set_dh_max(0.06)
        mdl.box_set_isotropic(0)
        mdl.alkane_set_dv_max(0.03  )



    else:
        mdl.alkane_set_dr_max(0.65)
        mdl.alkane_set_dt_max(0.43)
        mdl.alkane_set_dh_max(0.4)
        mdl.box_set_isotropic(1)
        mdl.alkane_set_dv_max(0.5)

    


   # Move types


    Sweeps_per_walk = 3*Nchains+1


    move_types = ['box','translate', 'rotate', 'dihedral']
    ivol = 0; itrans = 1; irot = 2; idih = 3
    moves_attempted = np.zeros(4)
    moves_accepted  = np.zeros(4)


    #adjusting mc step
    mc_adjust_interval = 1000
    interval_moves_attempted = np.zeros(4)
    interval_moves_accepted  = np.zeros(4)

    # Use the one and only simulation box
    ibox = 1

    # Initialise counters
    isweep = 0
    isweep_equil = 0

    #Equilibration parameters
    equil_sweeps = 50000 #steps to do until fully equilibrated
    equil_moves_attempted = np.zeros(4)
    equil_moves_accepted  = np.zeros(4)
    sample_sweeps= 100000 #number of steps to do once equilibrated
    equil_check = 20000
    equil = False #flag for equilibration. If true, equilibration sweeps are finished
    equil_error = 25
    
    ###unautomate equilibration
    ##mc step adjustment size
    equil_factor = 1.2 #factor to multiply the the step size by
    dv_min = 1e-5 #minimum step sizes for various moves, during the MC adjustment
    dr_min = 1e-5
    dt_min = 1e-5
    dh_min = 1e-5
    

    Max_sweeps = equil_sweeps + sample_sweeps #total number of sweeps

    # How often to sample quantities
    sample_interval = 10
    samples = list()

    # In case we're using Verlet lists rather than link cells
    mdl.alkane_construct_neighbour_list(ibox)

    # Initialise timer
    t1 = mdl.timer_init()

    #Visualisation parameters
    vis_interval = 500

    print(f"STARTING LOOP AT {pressure}")




    # Loop over the number of MC sweeps to perform
    while (isweep < Max_sweeps):

        if isweep%vis_interval == 0:
            #io.write(f'{workingdir}/configs.{isweep//vis_interval}.extxyz', mk_ase_config(ibox,Nbeads,Nchains))
            traj.write(atoms = mk_ase_config(ibox,Nbeads,Nchains))
            #output extxyz as well


        if isweep == equil_sweeps:
            equil = True
            logfile.write(f"{isweep} equilibration steps have been ran. Sampling configurations. \n")
            if ((np.std(samples[-(equil_check//sample_interval):])) > equil_error):
                logfile.write(f"Warning, system may not be fully equilibrated after {isweep} steps. \n")
                
        # One "sweep" is usually interpretted as one attempt to change each degree of freedom on average once.
        # here we have 3 translation + 2 rotation + 1 internal degrees of freedom per chain, plus 6 degrees
        # of freedom for the simulation cell (ignoring rotations).
        
        # A box move changes 1 degree of freedom
        # A chain translation move changes 3 degrees of freedom
        # A chain rotation move changes 2 degrees of freedom
        # A dihedral angle move changes 1 degree of freedom
       
        # Hence if we attempt each move type with a probability proportional to the number of degrees
        # of freedom it changes, we need to do 2*Nchains+6 moves to get one "sweep". Sampling in this
        # ratio isn't necessarily optimal, but it's a good starting point.
        
        # Calculate cumulative move ratios used to decide move type
        total_deg = 6*Nchains+6

        box_move_prob   = 6/total_deg
        trans_move_prob = box_move_prob   + 3.0*Nchains/total_deg
        rot_move_prob   = trans_move_prob + 2.0*Nchains/total_deg
        dih_move_prob   = rot_move_prob   + 1.0*Nchains/total_deg


        if isweep % 10000 == 0:
            logfile.write(f"{isweep} steps have been performed for {pressure} pressure in state {state}\n")
            
            
        


         



        # Loop over move attempts at the current sweep
        imove = 0
        while imove < 2*Nchains+6:

            #==========================#
            # Make a random trial move #
            #==========================#
            
            # Pick a random chain numbered from 0 to Nchains
            ichain = np.random.randint(0, high=Nchains)

            # Backup old chain positions. Note that whenever we call a routine inside the
            # hs_alkane library we need to add one to the bead/chain index to account for
            # Fortran indexing. 
            current_chain = mdl.alkane_get_chain(ichain+1, ibox)
            backup_chain = current_chain.copy() # Note copy, not equivalence

            # Pick a move type at random and call the appropriate function to make that 
            # trial move and return the corresponding ratio of Boltzmann factors.
            xi = np.random.random()
            if xi < box_move_prob:
                # Attempt a volume move
                itype = ivol
                boltz = mdl.alkane_box_resize(pressure, ibox, 0)
            elif xi < trans_move_prob:
                # Attempt a translation move
                itype = itrans
                boltz = mdl.alkane_translate_chain(ichain+1, ibox)
            elif xi < rot_move_prob:
                # Attempt a rotation move
                itype = irot
                boltz, quat = mdl.alkane_rotate_chain(ichain+1, ibox, 0)
            else:
                # Attempt a dihedral angle move
                itype = idih
                boltz, bead1, angle = mdl.alkane_bond_rotate(ichain+1, ibox, 1)
            
            # Increment attempted move counter
            moves_attempted[itype] += 1
            interval_moves_attempted[itype] +=1
            if equil:
                    equil_moves_attempted[itype] += 1
            
            #====================#
            # Accept/reject move #
            #====================#
            
            # Reject according to Metropolis criterion
            if (np.random.random() < boltz):
                
                # Move accepted
                moves_accepted[itype] += 1
                interval_moves_accepted[itype] += 1
                if equil:
                    equil_moves_accepted[itype] += 1
                
                # Update linked-list for new positions if not volume move. 
                if (itype!=ivol):
                    new_chain = mdl.alkane_get_chain(ichain+1, ibox).copy()
                    for ibead in range(Nbeads):
                        mdl.alkane_update_linked_lists(ibead+1, ichain+1, ibox, backup_chain[ibead], new_chain[ibead])
                        
            else:
                
                # Reject move
                if (itype!=ivol):
                    # Restore old chain if single chain move
                    for ibead in range(Nbeads):
                        current_chain[ibead] = backup_chain[ibead]
                else:
                    # Reset the box change - special fucntion for this.
                    dumboltz = mdl.alkane_box_resize(pressure, ibox, 1)

            imove += 1
            

            
        
            
        # Sample 
        if isweep%sample_interval == 0:
            samples.append(mdl.box_compute_volume(ibox))
            
#put a minimum step size

        if not equil:      #previously 'if True:'
        
            if (isweep%mc_adjust_interval == 0) and (isweep != 0) :
                equilfile.write(f"\nThis is step adjustment : {isweep//mc_adjust_interval}\n")
                if (interval_moves_accepted[0]/interval_moves_attempted[0]) > 0.50:
                    mdl.alkane_set_dv_max(mdl.alkane_get_dv_max()*equil_factor)
                elif (interval_moves_accepted[0]/interval_moves_attempted[0]) < 0.30:
                    mdl.alkane_set_dv_max(max(mdl.alkane_get_dv_max()/equil_factor,dv_min))
                else:
                    pass
                equilfile.write(f"Volume move acceptance rate is {interval_moves_accepted[0]/interval_moves_attempted[0]}\
                    {interval_moves_accepted[0]},{interval_moves_attempted[0]} \nVolume move size is now {mdl.alkane_get_dv_max()} \n")

                if (interval_moves_accepted[1]/interval_moves_attempted[1]) > 0.50:
                    mdl.alkane_set_dr_max(mdl.alkane_get_dr_max()*equil_factor)
                elif (interval_moves_accepted[1]/interval_moves_attempted[1]) < 0.30:
                    mdl.alkane_set_dr_max(max(mdl.alkane_get_dr_max()/equil_factor,dr_min))
                else:
                    pass
                equilfile.write(f"Translation move acceptance rate is {interval_moves_accepted[1]/interval_moves_attempted[1]}\
                    {interval_moves_accepted[1]},{interval_moves_attempted[1]} \nTranslation move size is now {mdl.alkane_get_dr_max()} \n")



                if (interval_moves_accepted[2]/interval_moves_attempted[2]) > 0.50:
                    mdl.alkane_set_dt_max(mdl.alkane_get_dt_max()*equil_factor)
                elif (interval_moves_accepted[2]/interval_moves_attempted[2]) < 0.30:
                    mdl.alkane_set_dt_max(max(mdl.alkane_get_dt_max()/equil_factor,dt_min))
                else:
                    pass
                equilfile.write(f"Angle move acceptance rate is {interval_moves_accepted[2]/interval_moves_attempted[2]} \
                    {interval_moves_accepted[2]},{interval_moves_attempted[2]} \nAngle move size is now {mdl.alkane_get_dt_max()} \n")


                if (interval_moves_accepted[3]/interval_moves_attempted[3]) > 0.50:
                    mdl.alkane_set_dh_max(mdl.alkane_get_dh_max()*equil_factor)
                elif (interval_moves_accepted[3]/interval_moves_attempted[3]) < 0.30:
                    mdl.alkane_set_dh_max(max(mdl.alkane_get_dh_max()/equil_factor,dh_min))
                else:
                    pass
                equilfile.write(f"Dihedral move acceptance rate is {interval_moves_accepted[3]/interval_moves_attempted[3]}\
                    {interval_moves_accepted[3]},{interval_moves_attempted[3]}, \nDihedral move size is now {mdl.alkane_get_dh_max()} \n")

               
                interval_moves_accepted[:] = 0
                interval_moves_attempted[:] = 0
        
        # Increment sweep counter and progress bar
        isweep += 1
        if equil:
            isweep_equil +=1
        
        # if isweep%sample_interval == 0:
        #     logfile.write(f"{isweep} for {pressure} pressure in {state} state")
        
        

    ##loop ends here

    traj.close()

        


    # Timing
    logfile.write(f"Completed at {pressure} pressure for {state} state {isweep} sweeps in {mdl.timer_elapsed_time()} seconds. \n")

    # Report statistics
    if equil:
        for itype in range(4):
            percent = 100.0 * equil_moves_accepted[itype]/equil_moves_attempted[itype]
            #print(itype, equil_moves_accepted[itype],equil_moves_attempted[itype], percent)
            logfile.write(f"Accepted {percent} % of {move_types[itype]} moves once equilibrated \n")
        
    # Rerun the sanity check to make sure nothing went wrong
    overlap_flag = mdl.alkane_check_chain_overlap(ibox)
    if overlap_flag != 0:
        logfile.write(f"Overlaps found between chains in configuration at pressure \n")
    

    else:
        logfile.write("No overlaps between chains found in configuration \n")
        
    equil_samples = samples[-(sample_sweeps//sample_interval):]


    # nlags = 200
    # x = np.array(equil_samples)
    # x = x - x.mean()
    # n = len(equil_samples)
    # autocovariance = np.empty(nlags + 1)
    # autocovariance[0] = x.dot(x)
    # for i in range(nlags):
    #    autocovariance[i + 1] = x[i + 1 :].dot(x[: -(i + 1)])
    # autocovariance /= n - np.arange(nlags + 1)
    # acf = autocovariance[: nlags + 1] /autocovariance[0]


    from statsmodels.tsa.stattools import acf
    Nlags = 2000
    simple_acf=acf(equil_samples, nlags=Nlags-1, fft=False)




    thinning_par = np.where(simple_acf==min(simple_acf))[0][0]
    
    logfile.write(f"Thinning parameter for decorrelation is {thinning_par} \n")
    logfile.write(f"Decorrelation value of {min(simple_acf)} \n")
    logfile.write(f"Final volume = {mdl.box_compute_volume(ibox)}\n")

    
    equilfile.close()
    equil_subsamples = equil_samples[::thinning_par]
    subsamples = samples[::thinning_par]

    np.savetxt(f"{workingdir}/volume.txt", subsamples)
    np.savetxt(f"{workingdir}/all_volumes.txt", samples)


    np.savetxt(f"{workingdir}/acf.txt", simple_acf)

    out = open(f"{maindir}/PvV_{state}.txt","a+")

    mean_vol = np.mean(equil_subsamples)
    stderr_vol = np.std(equil_subsamples)/np.sqrt(len(equil_subsamples))
    out.write("{},{},{}\n".format(pressure, mean_vol, stderr_vol))
    out.close()
    


    # Save volume samples
    
    logfile.write(f"Attempting to destroy simulation box \n")

    mdl.alkane_destroy()
    mdl.box_destroy_link_cells()
    mdl.box_destroy()

    logfile.write(f"Simulation box destroyed \n")

    logfile.close()    

###########################################
    
    
maindir = "./"

if not os.path.exists(maindir):
    os.makedirs(maindir)
with joblib.parallel_backend("multiprocessing"):
    Parallel(n_jobs=10)(delayed(alkane_MC)(pressure = i, solid = False, maindir = maindir) for i in range(10,50,5))

