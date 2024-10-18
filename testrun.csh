#! /bin/csh -f

./compile_all_mac.csh

rm -rf results_ant_template_test_1

mpiexec  -n 5  UFEMISM_program  config-files/config_ant_template_test.cfg