![Coverpage](https://github.com/wangshaoyun/Bond_Fluctuation_Method_SPME/blob/master/coverpage.jpg "Simulation System")
# Molecular Dynamics Simulation of Polyelectrolyte Brushes.
## About this program
1. This is a program used to simulate weak polyelectrolytes(PE) brushes with salts by Langevin dynamics [1,2].
2. This is the source program of Ref. [2-5].
3. Linear, star or mixing of linear and star PE brushes can be simulated. 
4. 8*8 chains are grafted on a plate, and each chain includes 100 monomers. Same quantities counterions of PE brushes are added in the system. In addition, trivalent salts with same charge quantities as charged monomers and corresponding couterions of trivalent salts are added in the system. The total charged particles are 21333. Moreover, external electric field can be considered.
5. Interaction includes Coulomb potential, Lenneard-Jones potential, finite extensiable nonlinear elastic (FENE) potential.
6. For short-range potential, the combination of cell lists method and Verlet lists method is used [6]. For long-range potential, smooth particle mesh Ewald (SPME) is used [7]. 
7. Slab geometry is used to correct the Ewald summation [8].
>[1] Ding H, Duan C, Tong C. Simulations of polymer brushes with charged end monomers under external electric fields, *Journal of Chemical Physics*, **146** (3), 034901, 2017.  
>[2] Zhang F, Wang S, Ding H, Tong C. Simulations of 3-arm polyelectrolyte star brushes under external electric fields. *Soft Matter*, **15** (12), 2560-2570, 2019.  
>[3] S. Y. Wang, C. H. Tong. Surface Switching of Mixed Polyelectrolyte Brushes Made of 4-arm Stars and Linear Chains: MD Simulations, *Journal of Applied Physics*, under review.  
>[4] T. B. Wang, S. Y. Wang, C. H. Tong. Charge Reversal of Polyelectrolyte Brushes Under a Collapsing Electric Field, to be submitted.  
>[5] Y. Ji, S. Y. Wang, C. H. Tong. The Collapse of Polyelectrolyte Brushes Made of 4-arm Stars Mediated by Trivalent Salt Ions and an Electric Field, to be submitted.  
>[6] Frenkel D, Klein M, Parrrinello M. Understanding Molecular Simulation: From Algorithm to Applications. 2nd ed. Singapore: Elsevier Pte Ltd, 2002. 
>[7] Ulrich Essmann, Lalith Perera, Max L. Berkowitz, Tom Darden, Hsing Lee, Lee G Pedersen. "A Smooth Particle Mesh Ewald Method." *Journal of Chemical Physics*, **103** (19), pp. 8577-8593, 1995.  
>[8] Yeh I C, Berkowitz M L. Ewald summation for systems with slab geometry. *Journal of Chemical Physics*, **111** (7), 3155-3162, 1999.

## About Parameter Files 
+ force_data.txt: It is a parameter file which includes Bjerrum length, external electric field
+ system_data.txt: It is a parameter file which includes system, time, brushes parameter.  

## Compiling files
1. Determine the parameters such as box size, brush chians, time step and so on.
2. Set the parameters in energy_data.txt and system_data.txt
3. Open terminal and enter the current file content
4. To compile the files:
```
$ make
```
  
5. To execute the files:
```
$ ./main &
```







