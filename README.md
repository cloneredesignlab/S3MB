**Stochastic State-Space Model of the Brain (S3MB)**. S3MB discretizes the brain into voxels, and the state of each voxel is defined by brain tissue stiffness, oxygen, glucose, vasculature, dead cells, migrating cells and proliferating cells of various ploidies, and treatment conditions such as radiotherapy and chemotherapy with temozolomide (TMZ). Well established Fokker-Planck partial differential equations govern the distribution of resources (oxygen, glucose) and cells across voxels:

<img width="904" alt="S3MB_Equations" src="https://github.com/cloneredesignlab/S3MB/assets/54926313/a8538526-5354-423c-afa6-efeaee3d88d2">

Interacting elements of the S3MB model are described in Figure 1. N indicates the total number of live tumor cells across all ploidy compartments n_i. Cell migration (K) is determined by brain tissue stiffness and by glucose and oxygen availability. In addition to influencing cell migration, oxygen availability has several downstream effects on the behaviour of the system: it impacts both chemo- and radiation therapy induced cell death (via parameters α_Tchemo and r* respectively) as well as switching between oxidative phosphorylation and glycolysis (further referred to as metabolic state). The metabolic state in turn also determines the rate at which cells proliferate (Michaelis-Menten constants α_g or β_g), die (Michaelis-Menten constants α_d or β_d) and consume glucose (rates ι_g and ι_g2). Both dead and live cell concentrations in focal and neighboring voxels determine the probability of vascular recruitment (p(V)). Stiffness of the brain tissue within each voxel (S) in conjunction with available resources dictates the migration potential of cells and is influenced by both treatment and tumor cell density.  

<img width="1281" alt="S3MB_ModelElements" src="https://github.com/cloneredesignlab/S3MB/assets/54926313/3dfe22f1-c68d-4daa-b411-e289929f4552">