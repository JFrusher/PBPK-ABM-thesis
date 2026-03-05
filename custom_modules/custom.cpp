/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#define _USE_MATH_DEFINES
#include <cmath>

#include "./custom.h"
#include <algorithm>   
#include <fstream>      
#include <sstream>      
#include <vector>
#include <iomanip>
#include <cstdlib>
#include <atomic>
#include <unordered_map>

namespace damage_params {
    double base_damage_rate = 0.001;          // Damage per uM drug per minute
    double s_phase_multiplier = 10.0;         // S-phase sensitivity multiplier
    double damage_arrest_threshold = 0.3;     // G1/S arrest if damage > 30%
    double damage_apoptosis_threshold = 0.7;  // Apoptosis if damage > 70%
    double transporter_upregulation = 0.0001; // ABC transporter induction rate
    double damage_repair = 0.00001;           // Baseline repair rate
    double quiescence_entry_rate = 0.0001;    // Stochastic quiescence entry
    double quiescence_stress_entry = 0.001;   // Damage-driven quiescence entry
    double quiescence_exit_rate = 0.00005;    // Exit rate when stress resolves
}

namespace drug_metabolism {
    // FdUMP dynamics (TS inhibition)
    double fdump_formation_rate = 0.015;        // 1/min (60 sec time constant)
    double fdump_decay_rate = 0.0033;           // 1/min (5 hour half-life, ~300 min)
    double fdump_TS_inhibition_potency = 5.0;   // Fold-reduction in dTMP synthesis per uM FdUMP
    
    // FdUTP dynamics (DNA incorporation)
    double fdutp_formation_rate = 0.020;        // 1/min (slightly faster than FdUMP)
    double fdutp_DNA_incorporation_rate_G1 = 0.001;    // Very low (no replication)
    double fdutp_DNA_incorporation_rate_S = 0.05;      // HIGH during S-phase (50 min^-1)
    double fdutp_DNA_incorporation_rate_G2M = 0.002;   // Low (replication complete)
    double fdutp_decay_rate = 0.0033;           // Same 5-hr half-life as parent
    
    // FUTP dynamics (RNA incorporation)
    double futp_formation_rate = 0.010;         // Slower (depends on FPRT enzyme)
    double futp_RNA_damage_rate = 0.008;        // Continuous while present
    double futp_decay_rate = 0.005;             // 3.3 hour half-life
    
    // DNA damage mechanics
    double FdUMP_dUTP_damage_per_uM = 0.0005;   // dUTP incorporation damage (slow accumulation)
    double FdUTP_direct_damage_per_incorporation = 0.01;  // Direct strand breaks (rapid)
    double FUTP_RNA_damage_per_uM = 0.0003;     // RNA damage slower to manifest
    
    // Repair processes
    double BER_repair_rate = 0.01;              // 1/min (100 min half-life for FdUTP damage)
    double HR_repair_rate = 0.002;              // Much slower, requires cell cycle restart
    double spontaneous_repair_rate = 0.00001;   // Background repair
    
    // Thymidylate synthase parameters
    double TS_baseline_expression = 1.0;        // Relative units
    double TS_upregulation_rate = 0.0001;       // Response to FdUMP burden
    double TS_max_expression = 3.0;             // Can increase 3-fold (resistance limit)
    double dTMP_synthesis_rate_normal = 50.0;   // uM/min (when TS active)
    
    // Checkpoint thresholds
    double DNA_damage_G1S_arrest_threshold = 0.2;      // 20% damage
    double DNA_damage_apoptosis_threshold = 0.6;       // 60% damage
    double DNA_damage_HR_activation_threshold = 0.4;   // Triggers HR (serious breaks)

}

// Multi-substrate PK profile loading
namespace pk_profiles {
    // Data containers for each substrate
    std::vector<double> time_points;           // Shared time axis
    
    // 5-FU concentration (from TUMOR_concentrations.csv)
    std::vector<double> concentration_5FU;
    
    // Metabolites (from METABOLITES.csv or PHYSICELL_METABOLITES.csv)
    std::vector<double> concentration_FdUMP;
    std::vector<double> concentration_FdUTP;
    std::vector<double> concentration_FUTP;
    std::vector<double> concentration_glucose;
    std::vector<double> concentration_thymidine;
    
    bool data_loaded = false;
    std::string loaded_prefix = "";
    
    // Mapping: which column in each CSV file contains substrate data
    // You'll set these based on your CSV column headers
    int column_5FU = -1;
    int column_folate = -1;
    int column_dNTP = -1;
    int column_glucose = -1;
    int column_thymidine = -1;
}

// Read a single column from a CSV file into a concentration series.
bool load_substrate_column(
    const std::string& filename,
    int column_index,
    int time_column,
    std::vector<double>& time_axis,
    std::vector<double>& concentration,
    bool time_is_minutes )
{
    std::ifstream file(filename);
    if( !file.is_open() )
    {
        std::cout << "WARNING: Could not open " << filename << std::endl;
        return false;
    }
    
    std::string line;
    std::vector<double> local_times;
    std::vector<double> local_conc;
    
    // Skip header
    std::getline(file, line);
    
    while( std::getline(file, line) )
    {
        std::stringstream ss(line);
        std::string token;
        std::vector<std::string> row;
        
        while( std::getline(ss, token, ',') )
        {
            row.push_back(token);
        }
        
        if( row.size() <= std::max(time_column, column_index) )
            continue;
        
        try {
            double t = std::stod(row[time_column]);
            double c = std::stod(row[column_index]);
            
            // Convert hours to minutes if needed
            if( !time_is_minutes )
                t *= 60.0;
            
            local_times.push_back(t);
            local_conc.push_back(c);
        }
        catch( ... ) {
            // Skip malformed rows
            continue;
        }
    }
    
    file.close();
    
    if( local_times.empty() )
        return false;
    
    // Initialize time axis from the first loaded substrate
    if( time_axis.empty() )
    {
        time_axis = local_times;
    }
    
    concentration = local_conc;
    return true;
}

// Load PK profiles for 5-FU and metabolites using a filename prefix.
void load_pk_profiles_multisubstrate(const std::string& prefix)
{
    std::cout << "\n========================================================" << std::endl;
    std::cout << "|    Loading 5-FU metabolite PK profiles                 |" << std::endl;
    std::cout << "========================================================" << std::endl;
    std::cout << "Prefix: " << prefix << std::endl;
    
    // Construct filenames for the selected dosing profile
    std::string file_tumor = prefix + "_TUMOR_concentrations.csv";
    std::string file_metabolites = prefix + "_PHYSICELL_METABOLITES.csv";
    
    std::cout << "\nLooking for files:" << std::endl;
    std::cout << "  1. " << file_tumor << std::endl;
    std::cout << "  2. " << file_metabolites << std::endl;
    
    // CSV column mapping
    // Column 0: Time_min
    // Column 2: Tumor_5FU_uM
    // Column 3: Tumor_FdUMP_uM
    // Column 4: Tumor_FdUTP_uM
    // Column 5: Tumor_FUTP_uM
    
    std::cout << "\n[1/4] Loading 5-FU (parent drug)..." << std::endl;
    if( load_substrate_column(file_tumor, 2, 1, pk_profiles::time_points, 
                               pk_profiles::concentration_5FU, false) )
    {
        std::cout << "      Loaded 5-FU: " << pk_profiles::concentration_5FU.size() << " points" << std::endl;
    }
    else
    {
        std::cout << "      Failed to load 5-FU" << std::endl;
    }
    
    std::cout << "[2/4] Loading FdUMP (TS inhibitor)..." << std::endl;
    if( load_substrate_column(file_metabolites, 3, 1, pk_profiles::time_points,
                               pk_profiles::concentration_FdUMP, false) )
    {
        std::cout << "      Loaded FdUMP: " << pk_profiles::concentration_FdUMP.size() << " points" << std::endl;
    }
    else
    {
        std::cout << "      Failed to load FdUMP" << std::endl;
    }
    
    std::cout << "[3/4] Loading FdUTP (DNA incorporation)..." << std::endl;
    if( load_substrate_column(file_metabolites, 4, 1, pk_profiles::time_points,
                               pk_profiles::concentration_FdUTP, false) )
    {
        std::cout << "      Loaded FdUTP: " << pk_profiles::concentration_FdUTP.size() << " points" << std::endl;
    }
    else
    {
        std::cout << "      Failed to load FdUTP" << std::endl;
    }
    
    std::cout << "[4/4] Loading FUTP (RNA incorporation)..." << std::endl;
    if( load_substrate_column(file_metabolites, 5, 1, pk_profiles::time_points,
                               pk_profiles::concentration_FUTP, false) )
    {
        std::cout << "      Loaded FUTP: " << pk_profiles::concentration_FUTP.size() << " points" << std::endl;
    }
    else
    {
        std::cout << "      Failed to load FUTP" << std::endl;
    }

    // Strict sanity check on required loaded vectors
    if( !pk_profiles::time_points.empty() )
    {
        if( pk_profiles::concentration_5FU.size() != pk_profiles::time_points.size() ||
            pk_profiles::concentration_FdUMP.size() != pk_profiles::time_points.size() ||
            pk_profiles::concentration_FdUTP.size() != pk_profiles::time_points.size() ||
            pk_profiles::concentration_FUTP.size() != pk_profiles::time_points.size() )
        {
            std::cout << "WARNING: PK profile length mismatch detected:" << std::endl;
            std::cout << "  time_points: " << pk_profiles::time_points.size() << std::endl;
            std::cout << "  5FU:         " << pk_profiles::concentration_5FU.size() << std::endl;
            std::cout << "  FdUMP:       " << pk_profiles::concentration_FdUMP.size() << std::endl;
            std::cout << "  FdUTP:       " << pk_profiles::concentration_FdUTP.size() << std::endl;
            std::cout << "  FUTP:        " << pk_profiles::concentration_FUTP.size() << std::endl;
        }
    }
    
    if( !pk_profiles::time_points.empty() )
    {
        pk_profiles::data_loaded = true;
        pk_profiles::loaded_prefix = prefix;
        
        std::cout << "\n========================================================" << std::endl;
        std::cout << "                    LOAD SUCCESSFUL                      " << std::endl;
        std::cout << "========================================================" << std::endl;
        std::cout << "\nPK Profile Summary:" << std::endl;
        std::cout << "  Time range: " << pk_profiles::time_points.front() << " - " 
                  << pk_profiles::time_points.back() << " min" << std::endl;
        std::cout << "  Data points: " << pk_profiles::time_points.size() << std::endl;
        std::cout << "\nSubstrate Status:" << std::endl;
        std::cout << "  5-FU:  LOADED" << std::endl;
        std::cout << "  FdUMP: LOADED" << std::endl;
        std::cout << "  FdUTP: LOADED" << std::endl;
        std::cout << "  FUTP:  LOADED" << std::endl;
    }
    else
    {
        std::cout << "\nWARNING: No data loaded." << std::endl;
    }
}


// Interpolate concentration by time for a given substrate series.
double get_concentration_at_time(const std::vector<double>& concentrations, 
                                  double current_time_min)
{
    if( pk_profiles::time_points.empty() || concentrations.empty() )
        return 0.0;
    
    if( current_time_min <= pk_profiles::time_points.front() )
        return concentrations.front();
    
    if( current_time_min >= pk_profiles::time_points.back() )
        return concentrations.back();
    
    // Linear interpolation
    for( size_t i = 0; i < pk_profiles::time_points.size() - 1; i++ )
    {
        if( current_time_min >= pk_profiles::time_points[i] && 
            current_time_min <= pk_profiles::time_points[i+1] )
        {
            double t0 = pk_profiles::time_points[i];
            double t1 = pk_profiles::time_points[i+1];
            double c0 = concentrations[i];
            double c1 = concentrations[i+1];
            
            double frac = (current_time_min - t0) / (t1 - t0);
            return c0 + frac * (c1 - c0);
        }
    }
    
    return concentrations.back();
}

// Update all microenvironment voxels from the loaded PK profiles.
void update_microenvironment_from_pk(double current_time)
{
    // Guard on availability
    if( !pk_profiles::data_loaded || pk_profiles::time_points.empty() )
        return;
    
    // Cache substrate indices
    static int idx_5FU = -1, idx_FdUMP = -1, idx_FdUTP = -1, 
               idx_FUTP = -1, idx_glucose = -1, idx_thymidine = -1;
    
    if( idx_5FU < 0 )
    {
        idx_5FU = microenvironment.find_density_index("5FU");
        idx_FdUMP = microenvironment.find_density_index("FdUMP");
        idx_FdUTP = microenvironment.find_density_index("FdUTP");
        idx_FUTP = microenvironment.find_density_index("FUTP");
        idx_glucose = microenvironment.find_density_index("glucose");
        idx_thymidine = microenvironment.find_density_index("thymidine");
    }
    
    // Pull concentrations for the current time
    double C_5FU = get_concentration_at_time(pk_profiles::concentration_5FU, current_time);
    double C_FdUMP = get_concentration_at_time(pk_profiles::concentration_FdUMP, current_time);
    double C_FdUTP = get_concentration_at_time(pk_profiles::concentration_FdUTP, current_time);
    double C_FUTP = get_concentration_at_time(pk_profiles::concentration_FUTP, current_time);
    double C_glucose = get_concentration_at_time(pk_profiles::concentration_glucose, current_time);
    double C_thymidine = get_concentration_at_time(pk_profiles::concentration_thymidine, current_time);
    
    // Update all voxels
    #pragma omp parallel for
    for( int i = 0; i < microenvironment.number_of_voxels(); i++ )
    {
        if( idx_5FU >= 0 )
            microenvironment(i)[idx_5FU] = C_5FU;
        if( idx_FdUMP >= 0 )
            microenvironment(i)[idx_FdUMP] = C_FdUMP;
        if( idx_FdUTP >= 0 )
            microenvironment(i)[idx_FdUTP] = C_FdUTP;
        if( idx_FUTP >= 0 )
            microenvironment(i)[idx_FUTP] = C_FUTP;
        // if( idx_glucose >= 0 )
        //     microenvironment(i)[idx_glucose] = C_glucose;
        // if( idx_thymidine >= 0 )
        //     microenvironment(i)[idx_thymidine] = C_thymidine;
    }
}

void load_pk_profile(std::string filename = "Chronomodulated_TUMOR_concentrations.csv");

// Initialize per-cell custom data used by the drug response model.
void init_damage_tracking(Cell* pCell)
{
    pCell->custom_data["FdUMP_burden"] = 0.0;           // uM*min of TS blockade
    pCell->custom_data["FdUTP_incorporated"] = 0.0;     // Cumulative FdUTP in DNA (0-1 scale)
    pCell->custom_data["FUTP_burden"] = 0.0;            // RNA damage accumulation
    
    // Phase and metabolic state tracking
    pCell->custom_data["current_phase"] = 0.0;          // 0=G1, 1=S, 2=G2, 3=M
    pCell->custom_data["metabolic_state"] = 0.5;        // 0=quiescent, 0.5=normal, 1=Warburg
    pCell->custom_data["folate_level"] = 100.0;         // Local availability
    pCell->custom_data["dNTP_level"] = 200.0;           // Local availability
    
    // Repair pathway activation
    pCell->custom_data["BER_activity"] = 0.0;           // Base excision repair (FdUTP damage)
    pCell->custom_data["HR_activity"] = 0.0;            // Homologous recombination (serious breaks)
    
    // Resistance markers
    pCell->custom_data["TS_expression"] = 1.0;          // Can be upregulated as resistance
    pCell->custom_data["transporter_MDR1"] = (pCell->type_name == "Cancer_stem") ? 0.5 : 0.0;
    pCell->custom_data["transporter_MRP"] = 0.0;        // MRP1 (exports metabolites)
    
    if (pCell->type_name == "Cancer_stem") {
        pCell->custom_data["ABC_transporter_expression"] = 0.5;
        pCell->custom_data["max_apoptosis_rate"] = 0.0001;
        pCell->custom_data["metabolic_state"] = 0.7;    // More quiescent baseline
    }
    if (pCell->type_name == "cancer_differentiated") {
        pCell->custom_data["max_apoptosis_rate"] = 0.005;
        pCell->custom_data["metabolic_state"] = 0.3;    // More Warburg-like
    }
}

// Configure cell definitions and bind custom callbacks.
void create_cell_types( void )
{
    // Set the random seed.
	if (parameters.ints.find_index("random_seed") != -1)
	{
		SeedRandom(parameters.ints("random_seed"));
	}
	
    // Initialize default cell definition and core callbacks.
	
	initialize_default_cell_definition(); 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	cell_defaults.functions.contact_function = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
    // Parse cell definitions from XML.
	initialize_cell_definitions_from_pugixml(); 

    // Build the map of cell definitions.
	build_cell_definitions_maps(); 

    // Initialize signal and response dictionaries.
	setup_signal_behavior_dictionaries(); 	

    // Load cell rules.
	setup_cell_rules(); 

    // Bind custom functions.
	cell_defaults.functions.update_phenotype = phenotype_function; 
	cell_defaults.functions.custom_cell_rule = custom_function; 
	cell_defaults.functions.contact_function = contact_function; 
	
    // Report setup.
	display_cell_definitions( std::cout ); 
	
	return; 
}

// Initialize the microenvironment.
void setup_microenvironment( void )
{
    // Initialize BioFVM.
	initialize_microenvironment(); 	
	
	return; 
}

// Randomize cell cycle phase by cell type.
void randomize_cell_phase_by_type( Cell* pC )
{
    int num_phases = 0;
    std::vector<double> phase_durations;
    
    if( pC->type_name == "Cancer_stem" )
    {
        num_phases = 4;
        phase_durations = {960.0, 480.0, 300.0, 60.0};  // G1, S, G2, M
    }
    else if( pC->type_name == "cancer_proliferating" )
    {
        num_phases = 4;
        phase_durations = {480.0, 360.0, 300.0, 60.0};  // G1, S, G2, M
    }
    else if( pC->type_name == "cancer_differentiated" )
    {
        num_phases = 2;
        phase_durations = {10000.0, 1000};  // Cycling quiescent
    }
    else
    {
        // Fallback
        num_phases = 3;
        phase_durations = {480.0, 360.0, 120.0};
    }
    
    // Randomize
    int random_phase = (int) floor( UniformRandom() * num_phases );
    pC->phenotype.cycle.data.current_phase_index = random_phase;
    
    double elapsed_time = UniformRandom() * phase_durations[random_phase];
    pC->phenotype.cycle.data.elapsed_time_in_phase = elapsed_time;
}

// Build initial tissue and configure oxygen Dirichlet nodes.
void setup_tissue( void )
{
    double Xmin = microenvironment.mesh.bounding_box[0]; 
    double Ymin = microenvironment.mesh.bounding_box[1]; 
    double Zmin = microenvironment.mesh.bounding_box[2]; 

    double Xmax = microenvironment.mesh.bounding_box[3]; 
    double Ymax = microenvironment.mesh.bounding_box[4]; 
    double Zmax = microenvironment.mesh.bounding_box[5]; 
    
    if( default_microenvironment_options.simulate_2D == true )
    {
        Zmin = 0.0; 
        Zmax = 0.0; 
    }
    
    double Xrange = Xmax - Xmin; 
    double Yrange = Ymax - Ymin; 
    double Zrange = Zmax - Zmin; 
    
    // Create initial cells for each type.
    
    Cell* pC;
    
    for( int k=0; k < cell_definitions_by_index.size() ; k++ )
    {
        Cell_Definition* pCD = cell_definitions_by_index[k]; 
        std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl; 
        for( int n = 0 ; n < parameters.ints("number_of_cells") ; n++ )
        {
            std::vector<double> position = {0,0,0}; 
            position[0] = Xmin + UniformRandom()*Xrange; 
            position[1] = Ymin + UniformRandom()*Yrange; 
            position[2] = Zmin + UniformRandom()*Zrange; 
            
            pC = create_cell( *pCD ); 
            pC->assign_position( position );
            init_damage_tracking(pC);
        }
    }
    std::cout << std::endl; 
    
    // Load cells from CSV and apply parameter distributions.
    load_cells_from_pugixml();
    set_parameters_from_distributions();
    
    // ========== HARD-STOP CAF DIVISION (ALTERNATIVE METHOD) ==========
    // Uncomment this block to completely disable CAF proliferation
    /*
    Cell_Definition* pCAF = find_cell_definition("CAF");
    if( pCAF != NULL )
    {
        std::cout << "BLOCKING CAF division: Setting all transition rates to 0.0" << std::endl;
        // Block all cycle transitions for CAF cells
        for( int i = 0; i < pCAF->phenotype.cycle.data.transition_rates.size(); i++ )
        {
            pCAF->phenotype.cycle.data.transition_rates[i] = 0.0;
        }
        // Alternatively, set a very long phase duration:
        // pCAF->phenotype.cycle.data.phase_durations[0] = 1e9;  // ~2 million years
    }
    */
    
    // Add Dirichlet nodes for oxygen supply.
    std::cout << "\n=== Setting up Dirichlet nodes for oxygen supply ===" << std::endl;

    int oxygen_index = microenvironment.find_density_index("oxygen");

    if( oxygen_index < 0 )
    {
        std::cout << "WARNING: Could not find oxygen substrate!" << std::endl;
    }
    else
    {
        std::cout << "Generating 'Core-Shell' vascular architecture for 3000 um tumor..." << std::endl;

        // Geometry
        double tumor_radius = 1600.0;        // 3000 um Diameter / 2
        double core_radius  = 400.0;         // Avascular / Necrotic Core
        double oxygen_value = 38.0;          // mmHg
        
        // Spacing parameters (denser at periphery)
        double icd_outer = 160.0;            // Intercapillary distance at margin
        double icd_inner = 300.0;            // Intercapillary distance near core
        double ring_spacing = 200.0;         // Distance between concentric shells
        
        int total_vessels = 0;
        srand(time(NULL)); 

        // Generate rings from core to tumor edge.
        for( double r = core_radius; r <= tumor_radius; r += ring_spacing )
        {
            // Linearly interpolate ICD based on radius.
            double frac = (r - core_radius) / (tumor_radius - core_radius);
            double current_icd = icd_inner - frac * (icd_inner - icd_outer);

            double circumference = 2.0 * M_PI * r;
            int num_vessels = (int)(circumference / current_icd);
            
            // Random phase shift so rings do not align.
            double phase_shift = (double)rand() / RAND_MAX * 2.0 * M_PI;

            std::cout << "  Ring at r=" << r << " um: Placing " << num_vessels 
                    << " vessels (ICD: " << (int)current_icd << " um)." << std::endl;

            for( int j = 0; j < num_vessels; j++ )
            {
                double angle = phase_shift + (2.0 * M_PI * j / num_vessels);
                
                // Jitter: +/- 30 um radial, +/- 0.05 radians angular.
                double r_jit = r + ((double)rand()/RAND_MAX * 60.0 - 30.0);
                double a_jit = angle + ((double)rand()/RAND_MAX * 0.1 - 0.05);

                // Keep jitter within bounds.
                if( r_jit < core_radius ) r_jit = core_radius;

                double x = r_jit * cos( a_jit );
                double y = r_jit * sin( a_jit );
                double z = 0.0; // 2D

                // Convert coordinates to voxel index.
                std::vector<double> vessel_position = {x, y, z};
                int voxel_idx = microenvironment.nearest_voxel_index(vessel_position);

                microenvironment.update_dirichlet_node(voxel_idx, oxygen_index, oxygen_value);
                total_vessels++;
            }
        }

        std::cout << "Successfully placed " << total_vessels << " vessels." << std::endl;
        std::cout << "  - Core Radius: " << core_radius << " um (Avascular)" << std::endl;
        std::cout << "  - Tumor Radius: " << tumor_radius << " um" << std::endl;
    }

    std::cout << "=== Dirichlet node setup complete ===" << std::endl;
    
    // Read dosing file prefix from environment variable, fallback to default.
    const char* dosing_file_env = std::getenv("DOSING_FILE");
    std::string dosing_file = (dosing_file_env != nullptr) ? dosing_file_env : "Chronomodulated";
    std::cout << "Using dosing file prefix: " << dosing_file << std::endl;
    load_pk_profiles_multisubstrate(dosing_file);
    
    return; 
}

// Coloring callback for SVG output.
std::vector<std::string> my_coloring_function( Cell* pCell )
{ return paint_by_number_cell_coloring(pCell); }

// Time-varying drug concentration from PK profile.

// PK profile storage.
std::vector<double> pk_time_points;      // Time in minutes
std::vector<double> pk_concentration;    // 5-FU concentration in uM
bool pk_data_loaded = false;

// Load a single-drug PK profile from CSV.
void load_pk_profile(std::string filename)
{
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cout << "ERROR: Could not open " << filename << std::endl;
        std::cout << "Using default constant concentration instead." << std::endl;
        return;
    }
    
    std::string line;
    std::getline(file, line);  // Skip header
    
    while (std::getline(file, line))
    {
        std::stringstream ss(line);
        std::string token;
        std::vector<std::string> row;
        
        while (std::getline(ss, token, ','))
        {
            row.push_back(token);
        }
        
        if (row.size() >= 3)
        {
            pk_time_points.push_back(std::stod(row[0]));
            pk_concentration.push_back(std::stod(row[2]));
        }
    }
    
    file.close();
    pk_data_loaded = true;
    
    std::cout << "========================================" << std::endl;
    std::cout << "PK PROFILE LOADED SUCCESSFULLY!" << std::endl;
    std::cout << "  Data points: " << pk_time_points.size() << std::endl;
    std::cout << "  Time range: 0 - " << pk_time_points.back() << " min" << std::endl;
    std::cout << "  Conc range: " << *std::min_element(pk_concentration.begin(), pk_concentration.end())
              << " - " << *std::max_element(pk_concentration.begin(), pk_concentration.end()) << " uM" << std::endl;
    std::cout << "========================================" << std::endl;
}

// Interpolate the single-drug PK profile.
double get_drug_concentration_at_time(double current_time_min)
{
    if (!pk_data_loaded || pk_time_points.empty())
    {
        return 0.0001;
    }
    
    // Linear interpolation between time points
    for (size_t i = 0; i < pk_time_points.size() - 1; i++)
    {
        if (current_time_min >= pk_time_points[i] && current_time_min < pk_time_points[i+1])
        {
            double t0 = pk_time_points[i];
            double t1 = pk_time_points[i+1];
            double c0 = pk_concentration[i];
            double c1 = pk_concentration[i+1];
            
            // Linear interpolation
            double fraction = (current_time_min - t0) / (t1 - t0);
            return c0 + fraction * (c1 - c0);
        }
    }
    
    return 0.0;
}
// Update microenvironment using the single-drug PK profile.
void update_5FU_in_microenvironment(double current_time)
{
    // Get current drug concentration from PK profile.
    double drug_conc = get_drug_concentration_at_time(current_time);
    
    // Get 5-FU substrate index.
    static int drug_index = microenvironment.find_density_index("5FU");
    
    // Update all voxels to have this concentration.
    #pragma omp parallel for
    for (int i = 0; i < microenvironment.number_of_voxels(); i++)
    {
        microenvironment(i)[drug_index] = drug_conc;
    }
    
//    static int apoptosis_sample_count = 0;
//	if (current_phase == 1 && apoptosis_sample_count < 20)  // Only S-phase cells
//	{
//		std::cout << "S-phase Cell " << pCell->ID 
//				<< " | Type: " << cell_type
//				<< " | Drug=" << drug_concentration 
//				<< " uM | Apoptosis=" << apoptosis_rate 
//				<< " /min (Death rate index [0]=" << pCell->phenotype.death.rates[0] << ")"
//				<< std::endl;
//		apoptosis_sample_count++;
//	}
    return;
}

// Main per-cell phenotype update (drug response and resistance).
void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{
    // int current_phase = pCell->phenotype.cycle.data.current_phase_index;
    // bool in_S_phase = (current_phase == PhysiCell_constants::Ki67_positive_premitotic);
    
    // // Adjust thymidine consumption based on phase
    // static int idx_thymidine = -1;
    // if( idx_thymidine == -1 )
    //     idx_thymidine = pCell->phenotype.secretion.find_secretion_substrate_index("thymidine");
    
    // if( idx_thymidine >= 0 )
    // {
    //     if( in_S_phase )
    //         pCell->phenotype.secretion.uptake_rates[idx_thymidine] = 0.005;  // Consume
    //     else
    //         pCell->phenotype.secretion.uptake_rates[idx_thymidine] = 0.0;     // Don't consume
    // }
    // Cache substrate indices.
    static int idx_5FU = microenvironment.find_density_index("5FU");
    static int idx_drug = microenvironment.find_density_index("drug");
    static int idx_FdUMP = microenvironment.find_density_index("FdUMP");
    static int idx_FdUTP = microenvironment.find_density_index("FdUTP");
    static int idx_FUTP = microenvironment.find_density_index("FUTP");

    // Abort if 5-FU is not found.
    if( idx_5FU < 0 ) return;

    // Read microenvironment concentrations.
    int voxelIndex = microenvironment.nearest_voxel_index( pCell->position );
    
    double C5FU = microenvironment.density_vector(voxelIndex)[idx_5FU];
    double C_drug = (idx_drug >= 0) ? microenvironment.density_vector(voxelIndex)[idx_drug] : 0.0;
    double C_FdUMP = (idx_FdUMP >= 0) ? microenvironment.density_vector(voxelIndex)[idx_FdUMP] : C5FU * 0.5;
    double C_FdUTP = (idx_FdUTP >= 0) ? microenvironment.density_vector(voxelIndex)[idx_FdUTP] : C5FU * 0.1;
    double C_FUTP  = (idx_FUTP >= 0) ? microenvironment.density_vector(voxelIndex)[idx_FUTP] : C5FU * 0.1;
    
    // Cell-cycle state.
    int current_phase = pCell->phenotype.cycle.data.current_phase_index; 
    // bool in_S_phase = (current_phase == PhysiCell_constants::Ki67_positive_premitotic);
    bool in_S_phase = (current_phase == 1);

    // Fetch resistance/response parameters (defaults if missing).
    auto get_custom_value = [&](const std::string& name, double fallback)
    {
        int idx = pCell->custom_data.find_variable_index(name);
        if( idx >= 0 )
            return pCell->custom_data[idx];
        return fallback;
    };

    double drug_sensitivity = get_custom_value("drug_sensitivity", 1.0);
    double abc_expression = get_custom_value("ABC_transporter_expression", 0.0);
    double efflux_5fu = get_custom_value("ABC_efflux_strength_5FU", 1.0);
    double efflux_met = get_custom_value("ABC_efflux_strength_metabolites", 2.0);

    if( abc_expression < 0.0 ) abc_expression = 0.0;
    if( abc_expression > 1.0 ) abc_expression = 1.0;

    double efflux_factor_5fu = 1.0 / (1.0 + efflux_5fu * abc_expression);
    double efflux_factor_met = 1.0 / (1.0 + efflux_met * abc_expression);

    double C5FU_eff = C5FU * efflux_factor_5fu;
    double C_FdUMP_eff = C_FdUMP * efflux_factor_met;
    double C_FdUTP_eff = C_FdUTP * efflux_factor_met;
    double C_FUTP_eff  = C_FUTP * efflux_factor_met;

    // Calculate damage rates.
    double FdUMP_damage_rate = 0.0;
    if( in_S_phase && C_FdUMP_eff > 0.0 )
    {
        FdUMP_damage_rate = drug_metabolism::FdUMP_dUTP_damage_per_uM * C_FdUMP_eff;
    }
    
    double FdUTP_damage_rate = (in_S_phase) ?
        (drug_metabolism::FdUTP_direct_damage_per_incorporation * C_FdUTP_eff) : 0.0;
    double FUTP_damage_rate  = drug_metabolism::FUTP_RNA_damage_per_uM * C_FUTP_eff;

    double total_damage_rate = (FdUMP_damage_rate + FdUTP_damage_rate + FUTP_damage_rate) * drug_sensitivity;

    // Apply damage to cell state.
    static int idx_damage = pCell->custom_data.find_variable_index("DNA_damage");
    if( idx_damage >= 0 )
    {
        double damage = pCell->custom_data[idx_damage];
        damage += total_damage_rate * dt;
        damage -= damage_params::damage_repair * damage * dt;
        if( damage < 0.0 ) damage = 0.0;
        if( damage > 1.0 ) damage = 1.0;
        pCell->custom_data[idx_damage] = damage;
    }

    // Apoptosis as a function of damage.
    double apoptosis_rate = 0.0;
    if( idx_damage >= 0 )
    {
        double damage = pCell->custom_data[idx_damage];
        double max_apoptosis_rate = get_custom_value("max_apoptosis_rate", 0.0);
        double threshold = damage_params::damage_apoptosis_threshold;
        if( damage > threshold )
        {
            double frac = (damage - threshold) / (1.0 - threshold);
            if( frac < 0.0 ) frac = 0.0;
            if( frac > 1.0 ) frac = 1.0;
            apoptosis_rate = max_apoptosis_rate * frac;
        }
    }

    int apoptosis_idx = phenotype.death.find_death_model_index( "apoptosis" );
    if( apoptosis_idx >= 0 )
    {
        phenotype.death.rates[apoptosis_idx] = apoptosis_rate;
    }

    // Resistance update (exposure-driven ABC induction).
    static int idx_exposure = pCell->custom_data.find_variable_index("cumulative_drug_exposure");
    if( idx_exposure >= 0 )
    {
        double exposure = pCell->custom_data[idx_exposure];
        exposure += (C5FU_eff + C_FdUMP_eff + C_FdUTP_eff + C_FUTP_eff) * dt;
        pCell->custom_data[idx_exposure] = exposure;
    }

    static int idx_abc = pCell->custom_data.find_variable_index("ABC_transporter_expression");
    if( idx_abc >= 0 )
    {
        double abc = pCell->custom_data[idx_abc];
        double induction = damage_params::transporter_upregulation *
            (C5FU_eff + C_FdUMP_eff + C_FdUTP_eff + C_FUTP_eff) * dt;
        abc += induction;
        if( abc < 0.0 ) abc = 0.0;
        if( abc > 1.0 ) abc = 1.0;
        pCell->custom_data[idx_abc] = abc;
    }

    // ========== CAF CONTACT INHIBITION ==========
    // Block CAF proliferation when local density is high
    if( pCell->type == 3 )  // CAF cell type ID
    {
        // Check neighbor density within ~50 microns (approximately 2-3 cell diameters)
        double interaction_radius = 50.0;  // Adjust this value to tune sensitivity
        int num_neighbors = pCell->cells_in_my_container().size() - 1;  // -1 to exclude self
        
        // Alternatively, use precise neighbor counting (slower but more accurate):
        // std::vector<Cell*> nearby = pCell->nearby_interacting_cells();
        // int num_neighbors = nearby.size();
        
        // Contact inhibition threshold: stop division if too crowded
        int neighbor_threshold = 8;  // Adjust this: higher = more permissive
        
        if( num_neighbors > neighbor_threshold )
        {
            // Find the G0/G1 -> S transition and block it
            // CAF uses cycling quiescent model (code 7) or flow cytometry (code 6)
            int cycle_code = pCell->phenotype.cycle.model().code;
            
            if( cycle_code == 7 )  // Cycling quiescent (Q -> G1 -> S)
            {
                // Block Q->G1 transition (index 0->1)
                pCell->phenotype.cycle.data.transition_rate(0,1) = 0.0;
            }
            else if( cycle_code == 6 )  // Flow cytometry (G0/G1 -> S -> G2 -> M)
            {
                // Block G0/G1->S transition (phase 0 -> phase 1)
                pCell->phenotype.cycle.data.transition_rate(0,1) = 0.0;
            }
            
            // Optional: Force cells already in cycle back to quiescence
            // (Uncomment if you want aggressive contact inhibition)
            // if( current_phase == 0 )  // If in G0/G1 phase
            // {
            //     pCell->phenotype.cycle.data.transition_rate(0,1) = 0.0;
            // }
        }
        else
        {
            // If not crowded, restore normal transition rate from XML settings
            // (PhysiCell automatically resets this from the cell definition XML each timestep,
            //  so you may not need to do anything here unless you're overriding frequently)
            // To explicitly restore: get the baseline rate from cell_definitions_by_type[3]
        }
    }

    return;
}

// Custom rule callback; used for diagnostics and PK updates.
void monitor_oxygenation( void )
{
    int oxygen_index = microenvironment.find_density_index("oxygen");
    int glucose_index = microenvironment.find_density_index("glucose");

    double current_time = PhysiCell_globals.current_time;

    double min_oxygen = 9e99;
    double max_oxygen = -9e99;
    double sum_oxygen = 0.0;

    double min_glucose = 9e99;
    double max_glucose = -9e99;
    double sum_glucose = 0.0;

    int n_vox = microenvironment.number_of_voxels();
    if( n_vox <= 0 )
    {
        std::cout << "[DIAG t=" << std::fixed << std::setprecision(1) << current_time
                  << " min] no voxels available" << std::endl;
        return;
    }

    for( int i = 0; i < n_vox; i++ )
    {
        if( oxygen_index >= 0 )
        {
            double oxygen_value = microenvironment(i)[oxygen_index];
            min_oxygen = std::min(min_oxygen, oxygen_value);
            max_oxygen = std::max(max_oxygen, oxygen_value);
            sum_oxygen += oxygen_value;
        }

        if( glucose_index >= 0 )
        {
            double glucose_value = microenvironment(i)[glucose_index];
            min_glucose = std::min(min_glucose, glucose_value);
            max_glucose = std::max(max_glucose, glucose_value);
            sum_glucose += glucose_value;
        }
    }

    double avg_oxygen = (oxygen_index >= 0) ? (sum_oxygen / n_vox) : 0.0;
    double avg_glucose = (glucose_index >= 0) ? (sum_glucose / n_vox) : 0.0;

    std::unordered_map<int, int> type_counts;
    int total_cells = 0;
    int alive_cells = 0;
    int dead_cells = 0;
    int apoptotic_cells = 0;
    int necrotic_cells = 0;

    for( Cell* pCell : *all_cells )
    {
        if( pCell == NULL )
        { continue; }

        total_cells++;
        type_counts[pCell->type]++;

        if( pCell->phenotype.death.dead )
        {
            dead_cells++;
            int death_model = pCell->phenotype.death.current_death_model_index;
            if( death_model == PhysiCell_constants::apoptosis_death_model )
            { apoptotic_cells++; }
            if( death_model == PhysiCell_constants::necrosis_death_model )
            { necrotic_cells++; }
        }
        else
        {
            alive_cells++;
        }
    }

    std::ostringstream type_summary;
    bool first_type = true;
    for( int i = 0; i < cell_definitions_by_type.size(); i++ )
    {
        auto* pDef = cell_definitions_by_type[i];
        if( pDef == NULL )
        { continue; }

        if( !first_type )
        { type_summary << " | "; }
        first_type = false;
        type_summary << pDef->name << "=" << type_counts[pDef->type];
    }

    if( first_type )
    { type_summary << "n/a"; }

    double C_5FU = get_concentration_at_time(pk_profiles::concentration_5FU, current_time);
    double C_FdUMP = get_concentration_at_time(pk_profiles::concentration_FdUMP, current_time);
    double C_FdUTP = get_concentration_at_time(pk_profiles::concentration_FdUTP, current_time);
    double C_FUTP = get_concentration_at_time(pk_profiles::concentration_FUTP, current_time);

    std::cout << "\n[DIAG t=" << std::fixed << std::setprecision(1) << current_time << " min]"
              << "------------------------------------------------------------" << std::endl;

    std::cout << "  ENV   | vox=" << n_vox;
    if( oxygen_index >= 0 )
    {
        std::cout << " | O2[min/avg/max]=" << std::setprecision(4)
                  << min_oxygen << "/" << avg_oxygen << "/" << max_oxygen;
    }
    else
    {
        std::cout << " | O2=n/a";
    }

    if( glucose_index >= 0 )
    {
        std::cout << " | Glucose[min/avg/max]="
                  << min_glucose << "/" << avg_glucose << "/" << max_glucose;
    }
    else
    {
        std::cout << " | Glucose=n/a";
    }
    std::cout << std::endl;

    std::cout << "  PK    | 5FU=" << C_5FU
              << " | FdUMP=" << C_FdUMP
              << " | FdUTP=" << C_FdUTP
              << " | FUTP=" << C_FUTP << std::endl;

    std::cout << "  CELLS | total=" << total_cells
              << " | alive=" << alive_cells
              << " | dead=" << dead_cells
              << " | apoptotic=" << apoptotic_cells
              << " | necrotic=" << necrotic_cells << std::endl;

    std::cout << "  TYPES | " << type_summary.str() << std::endl;
}

// Return whether the model should remain in diffusion-only warmup mode.
//
// Behavior:
// - During warmup, cell updates are paused in main.cpp (no cycle/death/mechanics/motility).
// - Diffusion-decay continues, allowing oxygen/glucose fields to equilibrate before biology starts.
//
// Configuration:
// - Optional XML user parameter: perfusion_warmup_minutes
// - If missing, default is 120 min.
// - Negative values are clamped to 0.
//
// Notes:
// - Parameter lookup is cached on first call.
// - A one-time console message is emitted when warmup is active.
bool perfusion_warmup_active( void )
{
    static bool initialized = false;
    static double warmup_minutes = 120.0; // Default: 2 hours diffusion-only equilibration.
    static bool reported = false;

    if( !initialized )
    {
        int idx = parameters.doubles.find_index("perfusion_warmup_minutes");
        if( idx >= 0 )
        {
            warmup_minutes = parameters.doubles("perfusion_warmup_minutes");
        }
        if( warmup_minutes < 0.0 )
        {
            warmup_minutes = 0.0;
        }
        initialized = true;
    }

    // Warmup remains active strictly before the configured cutoff time.
    bool active = (PhysiCell_globals.current_time < warmup_minutes);
    if( active && !reported )
    {
        std::cout << "[WARMUP] Diffusion-only mode active for first "
                  << warmup_minutes << " min (cell updates paused)." << std::endl;
        reported = true;
    }
    return active;
}

void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{
    double current_time = PhysiCell_globals.current_time;

    // This callback is executed for many cells (and potentially multiple threads).
    // Gate PK updates so only one caller updates once per simulation timestamp.
    static std::atomic<long long> last_update_stamp_ms(-1);
    long long current_stamp_ms = static_cast<long long>(std::llround(current_time * 1000.0));
    if( last_update_stamp_ms.exchange(current_stamp_ms) != current_stamp_ms )
    {
        update_microenvironment_from_pk(current_time);
    }

    // Emit diagnostics every 60 simulated minutes.
    int current_block = static_cast<int>( current_time / 60.0 );
    static std::atomic<int> last_monitor_block(-1);
    int seen_block = last_monitor_block.load();
    if( current_block > seen_block &&
        last_monitor_block.compare_exchange_strong(seen_block, current_block) )
    {
        monitor_oxygenation();
    }
    return;
}

// Contact function (no-op placeholder).
void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return; } 