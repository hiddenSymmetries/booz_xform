# Defaults for specific known hosts (i.e. HPC computing systems) can be put here. 

# For now, the following 4 lines are not used for anything. They might be useful in the future for detecting specific hosts.
cmake_host_system_information(RESULT HOST QUERY HOSTNAME)
cmake_host_system_information(RESULT FQDN QUERY FQDN)
message("hostname is ${HOST}")
message("fully qualified domain name is ${FQDN}")


if(DEFINED ENV{GITHUB_ACTIONS})
  message("Detected host is Github Actions CI")
  # Some flags are added to mpiexec to avoid some warning messages.
  set(BOOZ_COMMAND_TO_SUBMIT_JOB "mpiexec -n NUM_PROCS --mca btl_base_warn_component_unused 0")
  #set(BOOZ_COMMAND_TO_SUBMIT_JOB "mpiexec -n NUM_PROCS --mca btl_base_warn_component_unused 0 --mca orte_base_help_aggregate 0")

elseif(DEFINED ENV{NERSC_HOST})
  message("Detected host is NERSC Cori")
  set(CMAKE_CXX_COMPILER "CC")

elseif("$ENV{CLUSTER}" STREQUAL "DRACO")
  message("Detected host is IPP Draco")
  # CMake detects gnu compiler unless you specify Intel:
  set(CMAKE_CXX_COMPILER "mpiicpc")

elseif("$ENV{CLUSTER}" STREQUAL "COBRA")
  message("Detected host is IPP Cobra")
  # CMake detects gnu compiler unless you specify Intel:
  set(CMAKE_CXX_COMPILER "mpiicpc")

elseif("$ENV{CLUSTER}" STREQUAL "RAVEN")
  message("Detected host is IPP Raven")
  # CMake detects gnu compiler unless you specify Intel:
  set(CMAKE_CXX_COMPILER "mpiicpc")

else()
  message("This host is not a known one with specific rules.")

endif()
