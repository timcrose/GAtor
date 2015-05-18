[modules]
relaxation_module = amber_relaxation
initial_pool_relaxation_module = amber_relaxation
initial_pool_module = read_relax_user_structures
random_gen_module = random_crystal_1
comparison_module = mol_structure_comparison_3
selection_module = structure_selection_4
mutation_module = mutation_1
crossover_module = mol_crossover_4
initial_pool_comparison_module = mol_structure_comparison_3

[run_settings]
number_of_structures = 20
parallel_on_core = None ;  will run GA in parallel on one core
recover_from_crashes = False
verbose = True
number_of_replicas = 1
delta_convergence = .001

[control]
control_in_directory = control
initial_pool  = geometry_amber.key 
control_in_filelist = geometry_amber.key 

[random_gen]

[initial_pool]
number_of_processors = 1
user_structures_dir=/home/vama/soft/chem2/GAtor/src/test_integrate/100_random_sg
user2_structures_dir=/home/vama/soft/chem2/GAtor/src/test_integrate/10_random_sg

[comparison]
initial_vol = 155.0000
angle_up_bound = 120
angle_low_bound = 60
vol_decimal_tol = 0.3
fine_vol_tol = .01 

always_pass_initial_pool = False ; allows for bypassing initial pool test
energy_comparison = 1.0
dist_array_tolerance = 20 ; if structures are the same the dist sum will be 0...need to come back and fine-tune this 
energy_window = None ; may be None or int

[mutation]
number_of_permutations = 0
forbidden_elements = C H N O 

[TINKER]
path_to_tinker_executable = /home/vama/chem2/tinker/tinker/bin/xtalmin

[FHI-aims]
number_of_processors = 4 
path_to_aims_executable = /home/curtis/Codes/bin/aims.052014.mpi.x
initial_moment = hund

[crossover]
crossover_minimum_interface_distance = 1.0

[stoichiometry]
C = 4
H = 10
N = 2
O = 4

