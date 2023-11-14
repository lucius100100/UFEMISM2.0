#! /bin/csh -f

./compile_all_mac.csh

rm -rf results_ant_spinup_part1
mpiexec  -n 2  UFEMISM_program  config-files/config_antarctica_spinup_part1.cfg

#rm -rf results_ant_spinup_part2
#mpiexec  -n 2  UFEMISM_program  config-files/config_antarctica_spinup_part2.cfg

#rm -rf results_ant_retreat_DIVA
#mpiexec  -n 2  UFEMISM_program  config-files/config_antarctica_retreat_DIVA.cfg

#rm -rf results_ant_retreat_coupled
#mpiexec  -n 2  UFEMISM_program  config-files/config_antarctica_retreat_coupled.cfg

#rm -rf results_ant_retreat_BPA
#mpiexec  -n 2  UFEMISM_program  config-files/config_antarctica_retreat_BPA.cfg