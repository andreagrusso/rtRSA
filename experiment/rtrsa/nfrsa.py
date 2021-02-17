# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 09:18:25 2020

@author: Andrea Gerardo Russo, Biomedical Engineer
PhD candidate in Neuroscience
University of Salerno, Fisciano, Italy


rtRSA.
Class that creates a rtRSA object that can be used to perform real-time 
neurofeedback experiment using the real-time RSA techique (Russo et al., 
in preparation).

To create an rtRSA object it is needed to pass a name (str), 
the number of components (int) of the representational space (usually 2) and
the distance metric to use (str: pearson, euclidean, mean_diff)  

"""

import os
import numpy as np
from scipy.spatial.distance import euclidean
import numpy_indexed as npi
from sklearn.metrics.pairwise import euclidean_distances
import json



class rtRSA:

    """
    
    
    Define an rtRSA object by giving a name, the number of component of the 
    representational space and the distance metric to use for the estimation
    of the representational space.
    """    
    def __init__(self, name, n_comp, dist_metric):
        
        """
        Initialization of the class. The number of components
        of the RS and the distance metric are compulsory.
        
        name = base name for the ouptuts
        
        n_comp = Integer number higher than 2 (for a rtRSA-NF experiment
                 this number is 2)
        
                
        """
        
        
        self.n_comp = n_comp
        
        
        if dist_metric == 'pearson':
            self.dist_metric = 1
        if dist_metric == 'euclidean':
            self.dist_metric = 2
        if dist_metric == 'mean_diff':
            self.dist_metric = 3
           
        self.name = name
        
        #representational dissimilarity matrix
        self.RDM = []
        #coordinates of the base stimuli in the representational space (RS)
        self.RS_coords = []
        #matrix to project a new stimulus in the RS
        self.RS_inv_mat = []
        #t-values of the base stimuli (NrOfCoordinates X NrOfBaseStimuli)
        self.base_stimuli = []
        #coordinates of the ROI in the functional native space (FMR) from 
        #which the t-values of the base stimuli are extracted
        self.func_coords = []
        #name of the base stimuli
        self.conditions = []
        
        



    def __str__(self):
        print('rtRSA object.\n Name: ' + self.name + '\n')
        print('Nr of Representational Space components: '+ self.n_comp + '\n')
        print('Distance metric selected: ' + self.dist_metric)
    
    
    def load_RDM(self,rdm):
         
        """

        Parameters
        ----------
        rdm : numpy array
            N x N matrix that encodes the reciprocal dissimilarities among 
            the N base stimuli. The RDM represents the Representational Space 
            with its full dimensions.

        Returns
        -------
        None.

        """
        self.RDM = rdm
        
        

    def load_base_stimuli(self,base_stimuli):
         
        """

        Parameters
        ----------
        base_stimuli : Numpy matrix
            A matrix M x N with M number of the ROI voxels and N number of 
            base stimuli. This matrix contains the t-values relative to the
            base stimuli.

        Returns
        -------
        None.

        """
        self.base_stimuli = base_stimuli
        
        
        
    def load_inversion_matrix(self,matrix):
         
        """

        Parameters
        ----------
        matrix : Numpy matrix
            A matrix M x N with M that is the dimension of the RS (usually 2)
            and N number of base stimuli. This matrix contains a set of 
            values that can project a new set of t-values relative to a new
            stimulus into the RS loaded.

        Returns
        -------
        None.
        """

        self.RS_inv_mat = matrix
        
        if self.n_comp != min(matrix.shape):
            print('Please check the number of components selected and/or the inversion matrix')



    def load_RS_coords(self,coords):
         
        """
        Parameters
        ----------
        coords : Numpy matrix
            A N x M matrix with N number of base stimuli and M dimensions of 
            the RS (usually 2). These coordinates define the position of the
            base stimuli inside the RS.

        Returns
        -------
        None.

        """
        self.RS_coords = coords



    def load_func_coords(self,coords):
        
        """

        Parameters
        ----------
        coords : Numpy matrix
            3D coordinates of the voxels in the fucntional space. This set of
            coordinates are used to extract the t-values from each of the base
            stimuli of the localizer and it is used to as a reference for the
            functional coordinates of the VOI used in the NF session.

        Returns
        -------
        None.

        """
        self.func_coords = coords 
        
     
        
    def load_conditions(self,conditions):
        """
        

        Parameters
        ----------
        conditions : lsit
            List of the name of the base stimuli.

        Returns
        -------
        None.

        """
        self.conditions = conditions
        
     

    def saveAs(self,directory, basename):
        
        """
        
        Parameters
        ----------
        directory : string
            Ouptut directory.
        basename : string
            Basename for a .json file where there are stored the location
            of all the properties of the object.

        Returns
        -------
        None.


        
        """

        #save RDM
        np.savetxt(os.path.join(directory,basename+'.rdm'),self.RDM)

        #save RS coords
        np.savetxt(os.path.join(directory,basename+'.rsc'),self.RS_coords)
        
        #save t values base stimuli
        np.savetxt(os.path.join(directory,basename+'.tvals'),self.base_stimuli)
        
        #save inversion matrix
        np.savetxt(os.path.join(directory,basename+'.mat'),self.RS_inv_mat)        
        
        #save functional coords
        np.savetxt(os.path.join(directory,basename+'.fnc'),self.func_coords)        
        
        #save name of the base stimuli
        f=open(os.path.join(directory,basename+'.cnd'),'w')
        [f.write('%s\n' % item) for item in self.conditions]
        #create json file where the locations of all the rtRSA object properties
        #are stored
        data = {}
        
        
        data['metric'] = self.dist_metric
        data['n_comp'] = self.n_comp
        data['name'] = self.name
        data['rdm'] = os.path.join(directory,basename+'.rdm')
        data['base_stimuli'] = os.path.join(directory,basename+'.tvals')
        data['rs_coords'] = os.path.join(directory,basename+'.rsc')
        data['inv_mat'] = os.path.join(directory,basename+'.mat')
        data['func_coords'] = os.path.join(directory,basename+'.fnc')
        data['conditions'] = os.path.join(directory,basename+'.cnd')
        
        with open(os.path.join(directory,basename+'.json'), 'w') as outfile:
            json.dump(data, outfile,indent=4)
        


    def load(self,config_file):
        """
        This function load into a rtRSA object a set of existing data
        and rtRSA properties.

        Parameters
        ----------
        config_file : string
            Path of the .json file where are stored all the locations of the
            data of the created rtRSA.

        Returns
        -------
        None.

        """
        
        config = json.load(open(config_file))
        
       # input_name = input('Name of the input directory:\n')
        
        self.n_comp = config['n_comp']
        
        self.dist_metric = config['metric']
           
        self.name = config['name']
        
        #representational dissimilarity matrix
        self.RDM = np.loadtxt(config['rdm'])
        
        #coordinates of the base stimuli in the representational space (RS)
        self.RS_coords = np.loadtxt(config['rs_coords'])
        
        #matrix to project a new stimulus in the RS
        self.RS_inv_mat = np.loadtxt(config['inv_mat'])
        
        #t-values of the base stimuli (NrOfCoordinates X NrOfBaseStimuli)
        self.base_stimuli = np.loadtxt(config['base_stimuli'])
        
        #coordinates of the ROI in the functional native space (FMR) from 
        #which the t-values of the base stimuli are extracted
        self.func_coords = np.loadtxt(config['func_coords'])
        
        #import name of the base stimuli
        f = open(config['conditions'],'r')
        self.conditions = [line.rstrip() for line in f.readlines()]
                

    def cmdscale(self):
        
        """                                                                                       
        Classical multidimensional scaling (MDS)                                                  
                                                                                                   
        Parameters                                                                                
        ----------                                                                                
        self.RDM : (n, n) array                                                                          
            Symmetric distance matrix.
    
        self.n_comp : int
            Number of components requested
                                                                
                                                                                                   
        Returns                                                                                   
        -------                                                                                   
        self.RS_coords : (n, p) array                                                                          
            Configuration matrix. Each column represents a dimension. Only the                    
            p dimensions corresponding to positive eigenvalues of B are returned.                 
            Note that each dimension is only determined up to an overall sign,                    
            corresponding to a reflection.                                                        
                                                                                                   
        self.RS_inv_mat : (p,n) array                                                                            
            Matrix for the projection of a new stimulus.                                                                     
                                                                                                   
        """

        # Number of points                                                                        
        n = len(self.RDM)
     
        # Centering matrix                                                                        
        H = np.eye(n) - np.ones((n, n))/n
     
        # YY^T                                                                                    
        B = -H.dot(self.RDM**2).dot(H)/2
     
        # Diagonalize                                                                             
        evals, evecs = np.linalg.eigh(B)
     
        # Sort by eigenvalue in descending order                                                  
        idx   = np.argsort(evals)[::-1]
        evals = evals[idx]
        evecs = evecs[:,idx]
     
        #computing the matrix only for the positive-eigenvalued components
        w, = np.where(evals > 0)
        L  = np.diag(np.sqrt(evals[w]))
        V  = evecs[:,w]
        Y  = V.dot(L)
        
        
        
        #transformation matrix for only the first 2 components    
        self.RS_coords  = Y[:,:self.n_comp]
        self.RS_inv_mat  = evecs[:,:self.n_comp].T/(np.sqrt(evals[:self.n_comp]))[:,None]
        

        
    def createRS(self,data):
        
        """
        
        Return a Representation Dissimilarity Matrix estimated with a specific
        metric, the corresponding Representational Space (usuallly a 2D space) 
        a inversion mattrix that allows the projection of a new set of t-values
        in this Representational Space.
        
        Parameters
        ----------
        data : Numpy matrix
            This is a N x V matrix, with N number of base stimuli and V
            number of the functional voxels from which the t-values have been
            extracted..

        Returns
        -------
        RDM : Numpy matrix
            A matrix N x N with N that is the number of the base stimuli that
            encodes the Representational Space.
        This function return RDM, but it assigns also RDM to the corresponding
        internal attribute. 

        """

    
        data = data.T
        
    
        if self.dist_metric == 1:
            # Use correlation distance
            RDM = 1-np.corrcoef(data)
    
        elif self.dist_metric == 2:
            # Use Eucledian distance
            
            RDM = euclidean_distances(data,data,'euclidean')
    
        elif self.dist_metric == 3:
            # Use absolute activation difference
            means = np.mean(data,axis=1) # Determine mean activation per condition
            m, n = np.meshgrid(means,means) # Create all possible combinations
            RDM = abs(m-n) # Calculate difference for each combination
    
        self.RDM = RDM
                
        #calling the other internal function
        self.cmdscale()
        
        return RDM
        
        
    def target_positioning(self,new_stim):
    
        """
        This function estimate the distances between the t-values relative 
        to the base stimuli (estimated externally) and the t-values of the new
        stimulus. It returns the coordinates of the stimulus in the estimated RS 
        
        Parameters
        ----------
        new_stim : Numpy array
            vector that contains the t-values relative to the current
            stimulus (1xV, where V is the number of voxels)
            

        Returns
        -------
        coords[0],coords[1] : pair of float values
            x,y coordinates of the new stimulus in the embedded RS space
    
        """
        
        #the number of stimuli are the columns of the matrix
        NrOfBaseStim = np.shape(self.base_stimuli)[1]
    
        if self.dist_metric == 1:
            # Use the Pearson correlation and transform it in a distance
            distances = [1-np.corrcoef(self.base_stimuli[:,idx], new_stim)[0][1] 
            for idx in range(NrOfBaseStim)]
    
        if self.dist_metric == 2:
            # Use the classical euclidean distance
            distances = [euclidean(self.base_stimuli[:,idx],new_stim) 
            for idx in range(NrOfBaseStim)]
    
    
        if self.dist_metric == 3:
            # Use absolute activation difference
            distances = [abs(np.mean(self.base_stimuli[:,idx])-np.mean(new_stim)) 
            for idx in range(NrOfBaseStim)]
    

        #square the distances of the new stimuli
        sq_dist = np.square(distances)
    
        #compute mean of the square distances
        rdm_mean = np.sum(np.square(self.RDM),axis=1)/len(self.RDM)
        
        #compute the actual coords in the embedding space 
        coords = -0.5*(self.RS_inv_mat).dot((sq_dist - rdm_mean).T)
    
        return coords[0],coords[1],distances          



    def  match_coords(self,curr_coords):
        """
        This function takes in input the functional coordinates of 
        the current ROI used for the NF run. This set of coordinates
        would be different from the coordinates used to estimate
        the RDM as the subject could have performed the localizer
        task in a different scanning session. The aim of the function is 
        to erode/dilate the input set of coordinates to match its lenght
        with the coordinates used in the localizer

        Parameters
        ----------
        curr_coords : Numpy array
            3D coordinates in the native functional space of the current
            ROI.

        Returns
        -------
        nf_coords : Numpy array
            3D coordinates in the current native functional space that has
            the same lenght of the localizer coordinates and that need to 
            be used to extract the tvalues to apply the rtRSA.

        """

        #find centroid of the current coords
        c_trg = (np.sum(curr_coords,axis=0)/len(curr_coords)).reshape(1,-1) 
        
        #find the difference in the number of voxels
        delta = len(curr_coords) - len(self.func_coords)
        
        #find distances for the current coords from their centroid
        d_trg = euclidean_distances(curr_coords,c_trg)
        
        #if the the NF-ROI is bigger we need to remove voxels  
        if (delta > 0):
            
            print('THE CURRENT ROI IS BIGGER THAN THE ROI USED FOR THE RS!\n')
            print('THIS FUNCTION IS NOT FULLY TESTED!\n')
            print('STATISTICS COULD BE UNRELIABLE!\n')
            
            
            
            sort_idx = np.argsort(d_trg,axis=0)
            coords2del = sort_idx[-delta-1:-1]
            nf_coords = np.delete(curr_coords,coords2del,axis=0)
            nf_coords = nf_coords.astype(int)

            
            print('Most distant voxels deleted')
        
        #if the NF-ROI is smaller we need to add voxels    
        elif (delta < 0):
            
            print('THE CURRENT ROI IS SMALLER THAN THE ROI USED FOR THE RS!\n')
            print('THIS FUNCTION IS NOT FULLY TESTED!\n')
            print('STATISTICS COULD BE UNRELIABLE!\n')
            
            #find the min and the max of the coordinates over the coordinates
            bbox_edges = np.array([np.min(curr_coords,axis=0),
                             np.max(curr_coords,axis=0)])
            
            #create a bounding-box of coordinates
            bbox = np.mgrid[int(bbox_edges[0,0]):int(bbox_edges[1,0]+1),
                int(bbox_edges[0,1]):int(bbox_edges[1,1]+1),
                int(bbox_edges[0,2]):int(bbox_edges[1,2]+1)].reshape(3, -1).T
            
            #intersect the bounding-box and the existing coordinates
            surplus_voxels = npi.exclusive(bbox,
                                           curr_coords.astype(int),
                                           axis=0)
            
            #select among the possible coordinates the most distant from the centroid
            tmp_dist = euclidean_distances(surplus_voxels,c_trg)
            sort_idx = np.argsort(tmp_dist,axis=0)
            coords2add = sort_idx[:-delta]
            nf_coords = np.vstack((curr_coords.astype(int),
                                         surplus_voxels[coords2add,:].reshape(abs(delta),3)))
               
            print(delta,'voxels added to the ROI') 
        else:
            #the lucky situation when the number of coordinates are the same
            print('Same number of voxels! Lucky you!')
            nf_coords = curr_coords
        
        
        
        return nf_coords
