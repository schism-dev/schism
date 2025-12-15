!> A collection of standard varialbes used by QSim. Variable already definded 
!! in `module fabm_standard_variables` are excluded here.
module qsim_standard_variables
   use fabm_types 
   implicit none
   private
   
   public :: qsim_variables

   type :: type_qsim_variables
      ! =======================================================================
      ! surface variables
      ! =======================================================================
      type(type_horizontal_standard_variable) ::                                                                   &
         surface_longwave_flux_in_air = type_horizontal_standard_variable("surface_longwave_flux_in_air", "W/m2"), &
         surface_longwave_flux        = type_horizontal_standard_variable("surface_longwave_flux",        "W/m2"), &
         air_temperature              = type_horizontal_standard_variable("air_temperature",              "degC"), &
         water_level                  = type_horizontal_standard_variable("water_level",                  "m")

      ! =======================================================================
      ! bottom
      ! =======================================================================
      type(type_horizontal_standard_variable) ::                                                                                 &
         bottom_shortwave_flux               =type_horizontal_standard_variable("bottom_shortwave_flux",                "W/m2"), &
         bottom_shortwave_flux_in_sediment   =type_horizontal_standard_variable("bottom_shortwave_flux_in_sediment",    "W/m2"), &
         bottom_longwave_flux                =type_horizontal_standard_variable("bottom_longwave_flux",                 "W/m2"), &
         bottom_longwave_flux_in_sediment    =type_horizontal_standard_variable("bottom_longwave_flux_in_sediment",     "W/m2"), &
         bottom_photosynthetic_radiative_flux=type_horizontal_standard_variable("bottom_photosynthetic_radiative_flux", "W/m2"), &
         sediment_temperature                =type_horizontal_standard_variable("sediment_temperature",                 "degC")
      
      ! =======================================================================
      ! pelagic variables
      ! =======================================================================
      ! -----------------------------------------------------------------------
      ! algae
      ! -----------------------------------------------------------------------
      type(type_interior_standard_variable) ::                                                                                   &
         cyanobacteria_biomass                  = type_interior_standard_variable("cyanobacteria_biomass", "mg/l"),              &
         cyanobacteria_carbon                   = type_interior_standard_variable("cyanobacteria_carbon", "mgC/l"),              &
         cyanobacteria_carbon_content           = type_interior_standard_variable("cyanobacteria_carbon_content", "mgC/mg"),     &
         cyanobacteria_carbon_chlorophyll_ratio = type_interior_standard_variable("cyanobacteria_carbon_chlorophyll_ratio",      &
                                                                                  "mgC/mgChla"),                                 &
         cyanobacteria_chlorophyll_a            = type_interior_standard_variable("cyanobacteria_chlorophyll_a",                 & 
                                                                                  "mgChla/l"),                                   &
         cyanobacteria_nitrogen                 = type_interior_standard_variable("cyanobacteria_nitrogen",                      & 
                                                                                  "mgN/l"),                                      &
         cyanobacteria_nitrogen_content         = type_interior_standard_variable("cyanobacteria_nitrogen_content",              & 
                                                                                  "mgN/mg"),                                     &
         cyanobacteria_phosphor                 = type_interior_standard_variable("cyanobacteria_phosphor",                      & 
                                                                                  "mgP/l"),                                      &
         cyanobacteria_phosphor_content         = type_interior_standard_variable("cyanobacteria_phosphor_content",              & 
                                                                                  "mgP/mg")
                                   
      type(type_interior_standard_variable) ::                                                                                   &
         diatoms_biomass                   = type_interior_standard_variable("diatoms_biomass",                 "mg/l"      ),   &
         diatoms_carbon                    = type_interior_standard_variable("diatoms_carbon",                  "mgC/l"     ),   &
         diatoms_carbon_content            = type_interior_standard_variable("diatoms_carbon_content",          "mgC/mg"    ),   &
         diatoms_carbon_chlorophyll_ratio  = type_interior_standard_variable("diatoms_carbon_chloropyll_ratio", "mgC/mgChla"),   &
         diatoms_chlorophyll_a             = type_interior_standard_variable("diatoms_chlorophyll_a",           "mgChla/l"  ),   &
         diatoms_nitrogen                  = type_interior_standard_variable("diatoms_nitogen",                 "mgN/l"     ),   &   
         diatoms_nitrogen_content          = type_interior_standard_variable("diatoms_nitrogen_content",        "mgN/mg"    ),   &
         diatoms_phosphor                  = type_interior_standard_variable("diatoms_phosphor",                "mgP/l"     ),   &   
         diatoms_phosphor_content          = type_interior_standard_variable("diatoms_phosphor_content",        "mgP/mg"    ),   &
         diatoms_silicium                  = type_interior_standard_variable("diatoms_silicium",                "mgSi/l"    ),   &   
         diatoms_silicium_content          = type_interior_standard_variable("diatoms_silicium_content",        "mgSi/mg"   )
      
      type(type_interior_standard_variable) ::                                                                                   &
         green_algae_biomass                  = type_interior_standard_variable("green_alage_biomass",                  "mg/l"      ),  &
         green_algae_carbon                   = type_interior_standard_variable("green_algae_carbon",                   "mgC/l"     ),  &
         green_algae_carbon_content           = type_interior_standard_variable("green_algae_carbon",                   "mgC/mg"    ),  &
         green_algae_carbon_chlorophyll_ratio = type_interior_standard_variable("green_algae_carbon_chlorophyll_ratio", "mgC/mgChla"),  &
         green_algae_chlorophyll_a            = type_interior_standard_variable("green_algae_chlorophyll_a",            "mgChla/l"  ),  &
         green_algae_nitrogen                 = type_interior_standard_variable("green_algae_nitrogen",                 "mgN/mg"    ),  &
         green_algae_nitrogen_content         = type_interior_standard_variable("green_algae_nitrogen_content",         "mgN/l"     ),  &
         green_algae_phosphor                 = type_interior_standard_variable("green_algae_phosphor",                 "mgP/mg"    ),  &
         green_algae_phosphor_content         = type_interior_standard_variable("green_algae_phosphor_content",         "mgP/l"     )
      
      type(type_interior_standard_variable) ::                                                                              &
         algae_biomass     = type_interior_standard_variable("algae_biomass",     "mg/l",     aggregate_variable = .true.), &
         algae_chlorophyll = type_interior_standard_variable("algae_chlorophyll", "µgChla/l", aggregate_variable = .true.), &
         algae_carbon      = type_interior_standard_variable("algae_carbon",      "mgC/l",    aggregate_variable = .true.), &
         algae_nitrogen    = type_interior_standard_variable("algae_nitrogen",    "mgN/l",    aggregate_variable = .true.), &
         algae_phosphor    = type_interior_standard_variable("algae_phosphor",    "mgP/l",    aggregate_variable = .true.), &
         algae_silicium    = type_interior_standard_variable("algae_silicium",    "mgSi/l",   aggregate_variable = .true.)

      ! -----------------------------------------------------------------------
      ! detritus
      ! -----------------------------------------------------------------------
      type(type_interior_standard_variable) ::                                                                                            &
         readily_degradable_poc      = type_interior_standard_variable("readily_degradable_particulate_organic_carbon",         "mgC/l"), &
         hardly_degradable_poc       = type_interior_standard_variable("hardly_degradable_particulate_organic_carbon",          "mgC/l"), &
         readily_degradable_doc      = type_interior_standard_variable("readily_degradable_dissovled_organic_carbon_compounds", "mgC/l"), &
         hardly_degradable_doc       = type_interior_standard_variable("hardly_degradable_dissovled_organic_carbon_compounds",  "mgC/l"), &
         refractory_organic_carbon   = type_interior_standard_variable("refractory_organic_carbon_compounds",                   "mgC/l"), &
         monomeric_carbon            = type_interior_standard_variable("monomeric_carbon",                                      "mgC/l"), &
         detritus_carbon             = type_interior_standard_variable("detritus_carbon",                                       "mgC/l"), &
         detritus_nitrogen           = type_interior_standard_variable("detritus_nitrogen",                                     "mgN/l"), &
         detritus_phosphor           = type_interior_standard_variable("detritus_phosphor",                                     "mgP/l")                                                
      
      type(type_interior_standard_variable) ::  &
         suspended_particulate_matter = type_interior_standard_variable("suspended_particulate_matter",  "mg/l", aggregate_variable = .true.)

      ! ---------------------------------------------------------------------------------------------------------------------------
      ! zooplankton
      ! ---------------------------------------------------------------------------------------------------------------------------
      type(type_interior_standard_variable) ::                                                   &
         rotifer_biomass     = type_interior_standard_variable("rotifier_biomass",     "mg/l" ), &
         rotifer_individuals = type_interior_standard_variable("rotifier_individuals", "ind/l"), &
         rotifer_carbon      = type_interior_standard_variable("rotifier_carbon",      "mgC/l"), &
         rotifer_nitrogen    = type_interior_standard_variable("rotifier_nitrogen",    "mgN/l"), &
         rotifer_phosphor    = type_interior_standard_variable("rotifier_phosphor",    "mgP/l")
      
      type(type_interior_standard_variable) ::                                                                                        &
         zooplankton_biomass     = type_interior_standard_variable("zooplankton_biomass",     "mg/l",  aggregate_variable = .true.),  &
         zooplankton_individuals = type_interior_standard_variable("zooplankton_individuals", "ind/l", aggregate_variable = .true.),  &
         zooplankton_carbon      = type_interior_standard_variable("zooplankton_carbon",      "mgC/l", aggregate_variable = .true.),  &
         zooplankton_nitrogen    = type_interior_standard_variable("zooplankton_nitrogen",    "mgN/l", aggregate_variable = .true.),  &
         zooplankton_phosphor    = type_interior_standard_variable("zooplankton_phosphor",    "mgP/l", aggregate_variable = .true.)      
      
      ! ---------------------------------------------------------------------------------------------------------------------------
      ! 'nutrients'
      ! ---------------------------------------------------------------------------------------------------------------------------
      type(type_interior_standard_variable) ::                                                &
         ammonium          = type_interior_standard_variable("ammonium",           "mgN/l" ), &
         carbon_dioxid     = type_interior_standard_variable("carbon_dioxid",      "mgC/l" ), &
         dissoved_phosphor = type_interior_standard_variable("dissolved phosphor", "mgP/l" ), &
         nitrite           = type_interior_standard_variable("nitrite",            "mgN/l" ), &
         nitrate           = type_interior_standard_variable("nitrate",            "mgN/l" ), &
         oxygen            = type_interior_standard_variable("oxygen",             "mg/l"  ), &
         oxygen_saturation = type_interior_standard_variable("oxygen_saturation",  "%"     ), &
         silicium          = type_interior_standard_variable("silicium",           "mgSi/l")
      
      type(type_interior_standard_variable) ::                                                                                         &
         total_nitrogen           = type_interior_standard_variable("total_nitrogen",           "mgN/l", aggregate_variable = .true.), &
         total_organic_carbon     = type_interior_standard_variable("total_organic_carbon",     "mgC/l", aggregate_variable = .true.), &
         biological_oxygen_demand = type_interior_standard_variable("biological_oxygen_demand", "mg/l",  aggregate_variable = .true.), &
         chemical_oxygen_demand   = type_interior_standard_variable("chemical_oxygen_demand",   "mg/l",  aggregate_variable = .true.)
             
      ! ---------------------------------------------------------------------------------------------------------------------------
      ! bacteria
      ! ---------------------------------------------------------------------------------------------------------------------------
      type(type_interior_standard_variable) ::                                                         &
         heterotrophic_bacteria = type_interior_standard_variable("heterotrophic_bacteria", "mgC/l"),  &
         nitrobacter            = type_interior_standard_variable("nitrobacter",            "mg/l"),   &
         nitrosomonas           = type_interior_standard_variable("nitrosomonas",           "mg/l")
      
      ! ---------------------------------------------------------------------------------------------------------------------------
      ! radiation
      ! ---------------------------------------------------------------------------------------------------------------------------
      type(type_interior_standard_variable) ::                                                                                          &
         downwelling_longwave_flux = type_interior_standard_variable("downwelling_longwave_flux",  "W/m2"),                             &
         total_light_attenuation   = type_interior_standard_variable("total_light_attenuation",    "1/m", aggregate_variable = .true.)
          
      ! ---------------------------------------------------------------------------------------------------------------------------
      ! hydrodynamics
      ! ---------------------------------------------------------------------------------------------------------------------------
      type(type_interior_standard_variable) ::                                 &
         flow_velocity = type_interior_standard_variable("flow_velocity", "m/s")
      
      type(type_horizontal_standard_variable) ::                                                                  &
         friction_coefficient = type_horizontal_standard_variable("strickler_friction_coefficient", "m(1/3)/s") , &
         shear_velocity       = type_horizontal_standard_variable("shear_velocity",                 "m/s"),       &
         hydraulic_radius     = type_horizontal_standard_variable("hydraulic_radius",               "m"),         &
         total_depth          = type_horizontal_standard_variable("total_depth",                    "m")

      !---------------------------------------------------------------------------------------------------------------------------
      ! geometry
      ! ---------------------------------------------------------------------------------------------------------------------------
      type(type_horizontal_standard_variable) ::                               &
         bottom_area = type_horizontal_standard_variable("bottom_area", "m2"), &
         river_width = type_horizontal_standard_variable("river_width", "m"),  &
         zone_number = type_horizontal_standard_variable("zone_number", "-")
   end type type_qsim_variables

   type(type_qsim_variables), parameter :: qsim_variables = type_qsim_variables()

end module qsim_standard_variables
