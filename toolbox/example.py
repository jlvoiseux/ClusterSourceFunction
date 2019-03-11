import cluster
import point_charges
import comptonprofile

# Insert path of your Multiwfn executable
WFN_path = ""

# Create the cluster
clust = cluster.cluster(WFN_path)

# Add atoms to your cluster using an existing text file
clust.read_txt_atoms("C2H6.txt")

# Print the state of your cluster
print("La composition initiale du cluster est:")
print(clust)

# Build the input file used to solve the Hartree-Fock equations for the cluster using Gamess
clust.build_gamess_input("C2H6")

# Add interestpoints for the computation of the source function
clust.add_interestpoint(("C", 0), 1) # Atom
clust.add_interestpoint(("C",0), 2, ("C",1)) # Bond

# Compute the source contributions with quality "1"
clust.process_AIM_multiwfn("C2H6", ["1"])

# Compute the optimal charge of the point charge replacing ("C", 0)
pc = point_charges.point_charges(clust)
pc.optimal_ponctual_charge(("C",0))

# View the molecule in your browser
clust.visualize()
    
    
    
    

